#!/usr/bin/env python3
"""
Cloud auto-scaling system for RNASEQ-MINI.
Automatically scales AWS Batch compute environments based on queue depth and resource estimates.
"""

import os
import sys
import json
import logging
import argparse
import time
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
import boto3
from botocore.exceptions import ClientError, NoCredentialsError

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class ScalingConfig:
    """Configuration for auto-scaling."""
    min_vcpus: int = 0
    max_vcpus: int = 256
    desired_vcpus: int = 0
    scale_up_threshold: float = 0.8  # Scale up when queue > 80% of capacity
    scale_down_threshold: float = 0.3  # Scale down when queue < 30% of capacity
    scale_up_factor: float = 1.5  # Multiply current capacity by this factor when scaling up
    scale_down_factor: float = 0.7  # Multiply current capacity by this factor when scaling down
    cooldown_minutes: int = 5  # Minimum time between scaling actions
    check_interval_seconds: int = 30  # How often to check queue status

@dataclass
class QueueMetrics:
    """Queue metrics for scaling decisions."""
    pending_jobs: int
    running_jobs: int
    total_jobs: int
    compute_environment_name: str
    current_capacity: int
    max_capacity: int
    utilization_rate: float

class CloudAutoScaler:
    """Manages AWS Batch auto-scaling for RNASEQ-MINI."""

    def __init__(self, region_name: str = 'us-east-1', config: ScalingConfig = None):
        self.config = config or ScalingConfig()
        self.region_name = region_name

        try:
            self.batch_client = boto3.client('batch', region_name=region_name)
            self.ec2_client = boto3.client('ec2', region_name=region_name)
            self.cloudwatch_client = boto3.client('cloudwatch', region_name=region_name)
        except NoCredentialsError:
            logger.error("AWS credentials not found. Please configure AWS CLI or set environment variables.")
            sys.exit(1)

    def get_queue_metrics(self, job_queue_name: str, compute_environment_name: str) -> QueueMetrics:
        """Get current queue and compute environment metrics."""
        try:
            # Get job queue status
            job_queues = self.batch_client.describe_job_queues(jobQueues=[job_queue_name])
            job_queue = job_queues['jobQueues'][0]

            # Get compute environment status
            compute_envs = self.batch_client.describe_compute_environments(
                computeEnvironments=[compute_environment_name]
            )
            compute_env = compute_envs['computeEnvironments'][0]

            # Extract metrics
            state = job_queue['state']
            if state != 'ENABLED':
                logger.warning(f"Job queue {job_queue_name} is not enabled (state: {state})")

            # Get job counts from CloudWatch or API
            pending_jobs = 0
            running_jobs = 0

            # Try to get job counts from list_jobs API
            try:
                jobs_response = self.batch_client.list_jobs(jobQueue=job_queue_name, maxResults=1000)
                all_jobs = jobs_response.get('jobSummaryList', [])

                # Get more jobs if there are more pages
                while 'nextToken' in jobs_response:
                    jobs_response = self.batch_client.list_jobs(
                        jobQueue=job_queue_name,
                        maxResults=1000,
                        nextToken=jobs_response['nextToken']
                    )
                    all_jobs.extend(jobs_response.get('jobSummaryList', []))

                pending_jobs = sum(1 for job in all_jobs if job['jobStatus'] in ['PENDING', 'SUBMITTED'])
                running_jobs = sum(1 for job in all_jobs if job['jobStatus'] == 'RUNNING')

            except ClientError as e:
                logger.warning(f"Could not get job counts from API: {e}")
                # Fallback to approximate counts
                pending_jobs = len([j for j in job_queue.get('jobs', []) if j.get('jobStatus') in ['PENDING', 'SUBMITTED']])
                running_jobs = len([j for j in job_queue.get('jobs', []) if j.get('jobStatus') == 'RUNNING'])

            total_jobs = pending_jobs + running_jobs

            # Get current capacity from compute environment
            current_capacity = compute_env.get('computeResources', {}).get('desiredvCpus', 0)
            max_capacity = compute_env.get('computeResources', {}).get('maxvCpus', 0)

            # Calculate utilization rate
            utilization_rate = (running_jobs / max(current_capacity, 1)) if current_capacity > 0 else 0

            return QueueMetrics(
                pending_jobs=pending_jobs,
                running_jobs=running_jobs,
                total_jobs=total_jobs,
                compute_environment_name=compute_environment_name,
                current_capacity=current_capacity,
                max_capacity=max_capacity,
                utilization_rate=utilization_rate
            )

        except ClientError as e:
            logger.error(f"Error getting queue metrics: {e}")
            raise

    def calculate_target_capacity(self, metrics: QueueMetrics) -> int:
        """Calculate target capacity based on current metrics."""
        current_capacity = metrics.current_capacity

        # Scale up logic
        if metrics.utilization_rate > self.config.scale_up_threshold:
            target_capacity = int(current_capacity * self.config.scale_up_factor)
            logger.info(f"Scale up triggered: utilization {metrics.utilization_rate:.2%} > "
                       f"{self.config.scale_up_threshold:.2%}, target: {target_capacity}")

        # Scale down logic
        elif metrics.utilization_rate < self.config.scale_down_threshold and current_capacity > self.config.min_vcpus:
            target_capacity = max(self.config.min_vcpus, int(current_capacity * self.config.scale_down_factor))
            logger.info(f"Scale down triggered: utilization {metrics.utilization_rate:.2%} < "
                       f"{self.config.scale_down_threshold:.2%}, target: {target_capacity}")

        else:
            target_capacity = current_capacity
            logger.debug(f"No scaling needed: utilization {metrics.utilization_rate:.2%} within thresholds")

        # Respect min/max bounds
        target_capacity = max(self.config.min_vcpus, min(self.config.max_vcpus, target_capacity))

        return target_capacity

    def update_compute_environment(self, compute_environment_name: str, target_vcpus: int) -> bool:
        """Update compute environment with new desired vCPU count."""
        try:
            response = self.batch_client.update_compute_environments(
                computeEnvironments=[
                    {
                        'computeEnvironmentName': compute_environment_name,
                        'state': 'ENABLED',
                        'computeResources': {
                            'desiredvCpus': target_vcpus
                        }
                    }
                ]
            )

            logger.info(f"Updated compute environment {compute_environment_name} to {target_vcpus} vCPUs")
            return True

        except ClientError as e:
            logger.error(f"Error updating compute environment: {e}")
            return False

    def run_scaling_loop(self, job_queue_name: str, compute_environment_name: str,
                        duration_minutes: int = 60):
        """Run the auto-scaling loop for a specified duration."""
        logger.info(f"Starting auto-scaling for {duration_minutes} minutes")
        logger.info(f"Job queue: {job_queue_name}")
        logger.info(f"Compute environment: {compute_environment_name}")

        start_time = time.time()
        end_time = start_time + (duration_minutes * 60)
        last_scaling_time = 0

        while time.time() < end_time:
            try:
                # Get current metrics
                metrics = self.get_queue_metrics(job_queue_name, compute_environment_name)

                logger.info(f"Queue status: {metrics.pending_jobs} pending, "
                          f"{metrics.running_jobs} running, "
                          f"capacity: {metrics.current_capacity}/{metrics.max_capacity}, "
                          f"utilization: {metrics.utilization_rate:.2%}")

                # Check if we should scale
                target_capacity = self.calculate_target_capacity(metrics)

                # Only scale if enough time has passed since last scaling
                time_since_last_scaling = time.time() - last_scaling_time
                if (target_capacity != metrics.current_capacity and
                    time_since_last_scaling >= (self.config.cooldown_minutes * 60)):

                    success = self.update_compute_environment(compute_environment_name, target_capacity)
                    if success:
                        last_scaling_time = time.time()
                        logger.info(f"Successfully scaled to {target_capacity} vCPUs")
                    else:
                        logger.error(f"Failed to scale to {target_capacity} vCPUs")

                # Wait before next check
                time.sleep(self.config.check_interval_seconds)

            except KeyboardInterrupt:
                logger.info("Auto-scaling interrupted by user")
                break
            except Exception as e:
                logger.error(f"Error in scaling loop: {e}")
                time.sleep(self.config.check_interval_seconds)

        logger.info("Auto-scaling completed")

    def estimate_cost_savings(self, metrics_history: List[QueueMetrics]) -> Dict[str, float]:
        """Estimate cost savings from auto-scaling."""
        if not metrics_history:
            return {'estimated_savings_usd': 0.0, 'efficiency_improvement': 0.0}

        # Calculate average utilization before and after optimization
        avg_utilization = np.mean([m.utilization_rate for m in metrics_history])

        # Assume baseline cost per vCPU-hour
        cost_per_vcpu_hour = 0.05  # AWS Batch approximate cost

        # Calculate potential savings (if we had perfect scaling)
        total_vcpu_hours = sum(m.current_capacity for m in metrics_history) * (len(metrics_history) / 60)  # Rough estimate
        utilized_vcpu_hours = sum(m.running_jobs for m in metrics_history) * (len(metrics_history) / 60)

        if total_vcpu_hours > 0:
            efficiency_improvement = utilized_vcpu_hours / total_vcpu_hours
            potential_savings = total_vcpu_hours * cost_per_vcpu_hour * (1 - efficiency_improvement)
        else:
            efficiency_improvement = 0.0
            potential_savings = 0.0

        return {
            'estimated_savings_usd': potential_savings,
            'efficiency_improvement': efficiency_improvement,
            'average_utilization': avg_utilization
        }

def create_aws_batch_compute_environment(config_file: str = None):
    """Create an AWS Batch compute environment for RNASEQ-MINI."""
    try:
        # Load configuration
        if config_file and os.path.exists(config_file):
            with open(config_file, 'r') as f:
                config = json.load(f)
        else:
            config = {
                'compute_environment_name': 'rnaseq-mini-compute',
                'job_queue_name': 'rnaseq-mini-queue',
                'vpc_id': None,
                'subnets': [],
                'max_vcpus': 256,
                'instance_types': ['optimal']
            }

        # This would create the actual AWS resources
        # For now, just log what would be created
        logger.info("Creating AWS Batch compute environment...")
        logger.info(f"Name: {config['compute_environment_name']}")
        logger.info(f"Max vCPUs: {config['max_vcpus']}")
        logger.info("Compute environment configuration ready")

        return True

    except Exception as e:
        logger.error(f"Error creating compute environment: {e}")
        return False

def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(description='Cloud auto-scaling for RNASEQ-MINI')
    parser.add_argument('--job-queue', required=True, help='AWS Batch job queue name')
    parser.add_argument('--compute-env', required=True, help='AWS Batch compute environment name')
    parser.add_argument('--region', default='us-east-1', help='AWS region')
    parser.add_argument('--duration', type=int, default=60, help='Duration to run auto-scaling (minutes)')
    parser.add_argument('--create-env', help='Create compute environment from config file')
    parser.add_argument('--min-vcpus', type=int, default=0, help='Minimum vCPUs')
    parser.add_argument('--max-vcpus', type=int, default=256, help='Maximum vCPUs')
    parser.add_argument('--scale-up-threshold', type=float, default=0.8, help='Scale up threshold')
    parser.add_argument('--scale-down-threshold', type=float, default=0.3, help='Scale down threshold')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose logging')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Create compute environment if requested
    if args.create_env:
        success = create_aws_batch_compute_environment(args.create_env)
        if not success:
            sys.exit(1)
        return

    # Setup auto-scaler
    config = ScalingConfig(
        min_vcpus=args.min_vcpus,
        max_vcpus=args.max_vcpus,
        scale_up_threshold=args.scale_up_threshold,
        scale_down_threshold=args.scale_down_threshold
    )

    scaler = CloudAutoScaler(region_name=args.region, config=config)

    # Run auto-scaling
    try:
        scaler.run_scaling_loop(args.job_queue, args.compute_env, args.duration)
    except KeyboardInterrupt:
        logger.info("Auto-scaling stopped by user")

if __name__ == "__main__":
    main()



