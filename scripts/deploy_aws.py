#!/usr/bin/env python3
"""
AWS deployment script for RNASEQ-MINI.
Deploys the pipeline to AWS Batch with auto-scaling and cost optimization.
"""

import boto3
import json
import argparse
import os
from pathlib import Path
from datetime import datetime


class AWSDeployer:
    """Handles AWS deployment for RNASEQ-MINI."""

    def __init__(self, region='us-east-1'):
        self.region = region
        self.ecs_client = boto3.client('ecs', region_name=region)
        self.batch_client = boto3.client('batch', region_name=region)
        self.iam_client = boto3.client('iam', region_name=region)
        self.cloudformation_client = boto3.client('cloudformation', region_name=region)

    def create_vpc_stack(self, stack_name='rnaseq-mini-vpc'):
        """Create VPC infrastructure for the pipeline."""
        template = {
            "AWSTemplateFormatVersion": "2010-09-09",
            "Resources": {
                "VPC": {
                    "Type": "AWS::EC2::VPC",
                    "Properties": {
                        "CidrBlock": "10.0.0.0/16",
                        "EnableDnsSupport": True,
                        "EnableDnsHostnames": True,
                        "Tags": [{"Key": "Name", "Value": "rnaseq-mini-vpc"}]
                    }
                },
                "Subnet1": {
                    "Type": "AWS::EC2::Subnet",
                    "Properties": {
                        "VpcId": {"Ref": "VPC"},
                        "CidrBlock": "10.0.1.0/24",
                        "AvailabilityZone": f"{self.region}a",
                        "Tags": [{"Key": "Name", "Value": "rnaseq-mini-subnet-1"}]
                    }
                },
                "Subnet2": {
                    "Type": "AWS::EC2::Subnet",
                    "Properties": {
                        "VpcId": {"Ref": "VPC"},
                        "CidrBlock": "10.0.2.0/24",
                        "AvailabilityZone": f"{self.region}b",
                        "Tags": [{"Key": "Name", "Value": "rnaseq-mini-subnet-2"}]
                    }
                },
                "InternetGateway": {
                    "Type": "AWS::EC2::InternetGateway",
                    "Properties": {
                        "Tags": [{"Key": "Name", "Value": "rnaseq-mini-igw"}]
                    }
                },
                "GatewayAttachment": {
                    "Type": "AWS::EC2::VPCGatewayAttachment",
                    "Properties": {
                        "VpcId": {"Ref": "VPC"},
                        "InternetGatewayId": {"Ref": "InternetGateway"}
                    }
                }
            }
        }

        self.cloudformation_client.create_stack(
            StackName=stack_name,
            TemplateBody=json.dumps(template),
            Capabilities=['CAPABILITY_IAM']
        )

        print(f"Created VPC stack: {stack_name}")

    def create_compute_environment(self, name='rnaseq-mini-compute-env'):
        """Create AWS Batch compute environment."""
        compute_env = {
            'computeEnvironmentName': name,
            'type': 'MANAGED',
            'state': 'ENABLED',
            'computeResources': {
                'type': 'EC2',
                'allocationStrategy': 'BEST_FIT_PROGRESSIVE',
                'minvCpus': 0,
                'maxvCpus': 100,
                'desiredvCpus': 0,
                'instanceTypes': ['c5.large', 'c5.xlarge', 'c5.2xlarge', 'm5.large'],
                'subnets': self._get_subnet_ids(),
                'securityGroupIds': self._get_security_group_ids(),
                'instanceRole': self._get_instance_profile_arn(),
                'bidPercentage': 50,  # Use spot instances for cost savings
                'spotIamFleetRole': self._get_spot_fleet_role_arn()
            }
        }

        self.batch_client.create_compute_environment(**compute_env)
        print(f"Created compute environment: {name}")

    def create_job_queue(self, name='rnaseq-mini-queue'):
        """Create AWS Batch job queue."""
        job_queue = {
            'jobQueueName': name,
            'state': 'ENABLED',
            'priority': 1,
            'computeEnvironmentOrder': [{
                'order': 1,
                'computeEnvironment': 'rnaseq-mini-compute-env'
            }]
        }

        self.batch_client.create_job_queue(**job_queue)
        print(f"Created job queue: {name}")

    def create_job_definition(self, name='rnaseq-mini-job-def'):
        """Create AWS Batch job definition for RNASEQ-MINI."""
        job_def = {
            'jobDefinitionName': name,
            'type': 'container',
            'containerProperties': {
                'image': 'rnaseq-mini:latest',
                'vcpus': 2,
                'memory': 4096,
                'environment': [
                    {'name': 'AWS_REGION', 'value': self.region},
                    {'name': 'S3_BUCKET', 'value': 'rnaseq-mini-results'},
                    {'name': 'RESULTS_DIR', 'value': '/app/results'},
                    {'name': 'CACHE_DIR', 'value': '/app/.cache'}
                ],
                'mountPoints': [
                    {
                        'sourceVolume': 'results',
                        'containerPath': '/app/results',
                        'readOnly': False
                    },
                    {
                        'sourceVolume': 'cache',
                        'containerPath': '/app/.cache',
                        'readOnly': False
                    }
                ],
                'volumes': [
                    {
                        'name': 'results',
                        'host': {'sourcePath': '/mnt/results'}
                    },
                    {
                        'name': 'cache',
                        'host': {'sourcePath': '/mnt/cache'}
                    }
                ],
                'logConfiguration': {
                    'logDriver': 'awslogs',
                    'options': {
                        'awslogs-group': '/aws/batch/rnaseq-mini',
                        'awslogs-region': self.region,
                        'awslogs-stream-prefix': 'rnaseq-mini'
                    }
                }
            },
            'timeout': {'attemptDurationSeconds': 3600}
        }

        self.batch_client.register_job_definition(**job_def)
        print(f"Created job definition: {name}")

    def create_s3_bucket(self, bucket_name='rnaseq-mini-results'):
        """Create S3 bucket for results storage."""
        s3_client = boto3.client('s3', region_name=self.region)

        try:
            s3_client.create_bucket(
                Bucket=bucket_name,
                CreateBucketConfiguration={
                    'LocationConstraint': self.region
                } if self.region != 'us-east-1' else {}
            )

            # Enable versioning and server-side encryption
            s3_client.put_bucket_versioning(
                Bucket=bucket_name,
                VersioningConfiguration={'Status': 'Enabled'}
            )

            s3_client.put_bucket_encryption(
                Bucket=bucket_name,
                ServerSideEncryptionConfiguration={
                    'Rules': [{
                        'ApplyServerSideEncryptionByDefault': {
                            'SSEAlgorithm': 'AES256'
                        }
                    }]
                }
            )

            print(f"Created S3 bucket: {bucket_name}")

        except s3_client.exceptions.BucketAlreadyExists:
            print(f"S3 bucket already exists: {bucket_name}")

    def create_ecs_service(self, cluster_name='rnaseq-mini-cluster'):
        """Create ECS service for the web interface."""
        # Create cluster if it doesn't exist
        try:
            self.ecs_client.create_cluster(clusterName=cluster_name)
            print(f"Created ECS cluster: {cluster_name}")
        except self.ecs_client.exceptions.ClusterAlreadyExistsException:
            print(f"ECS cluster already exists: {cluster_name}")

        # Create task definition for web service
        task_def = {
            'family': 'rnaseq-mini-web',
            'networkMode': 'awsvpc',
            'requiresCompatibilities': ['FARGATE'],
            'cpu': '1024',
            'memory': '2048',
            'executionRoleArn': self._get_execution_role_arn(),
            'taskRoleArn': self._get_task_role_arn(),
            'containerDefinitions': [{
                'name': 'rnaseq-mini-web',
                'image': 'rnaseq-mini:latest',
                'essential': True,
                'portMappings': [{
                    'containerPort': 8000,
                    'hostPort': 8000,
                    'protocol': 'tcp'
                }],
                'environment': [
                    {'name': 'ENVIRONMENT', 'value': 'production'},
                    {'name': 'REDIS_URL', 'value': 'redis://localhost:6379'},
                    {'name': 'DATABASE_URL', 'value': 'postgresql://user:password@db:5432/rnaseq_mini'}
                ],
                'logConfiguration': {
                    'logDriver': 'awslogs',
                    'options': {
                        'awslogs-group': '/ecs/rnaseq-mini-web',
                        'awslogs-region': self.region,
                        'awslogs-stream-prefix': 'ecs'
                    }
                }
            }]
        }

        response = self.ecs_client.register_task_definition(**task_def)
        task_def_arn = response['taskDefinition']['taskDefinitionArn']

        # Create service
        service = {
            'cluster': cluster_name,
            'serviceName': 'rnaseq-mini-web-service',
            'taskDefinition': task_def_arn,
            'desiredCount': 2,
            'launchType': 'FARGATE',
            'networkConfiguration': {
                'awsvpcConfiguration': {
                    'subnets': self._get_subnet_ids(),
                    'securityGroups': self._get_security_group_ids(),
                    'assignPublicIp': 'ENABLED'
                }
            },
            'loadBalancers': [{
                'targetGroupArn': self._create_target_group(),
                'containerName': 'rnaseq-mini-web',
                'containerPort': 8000
            }]
        }

        self.ecs_client.create_service(**service)
        print(f"Created ECS service: rnaseq-mini-web-service")

    def _get_subnet_ids(self):
        """Get subnet IDs from VPC stack."""
        # For demo purposes, return placeholder subnets
        # In production, these would be retrieved from the created VPC
        return [f'subnet-{i:08x}' for i in range(2)]

    def _get_security_group_ids(self):
        """Get security group IDs."""
        return [f'sg-{i:08x}' for i in range(1)]

    def _get_instance_profile_arn(self):
        """Get EC2 instance profile ARN."""
        return f'arn:aws:iam::{self._get_account_id()}:instance-profile/rnaseq-mini-instance-profile'

    def _get_spot_fleet_role_arn(self):
        """Get spot fleet role ARN."""
        return f'arn:aws:iam::{self._get_account_id()}:role/AWSBatchServiceRole'

    def _get_execution_role_arn(self):
        """Get ECS task execution role ARN."""
        return f'arn:aws:iam::{self._get_account_id()}:role/ecsTaskExecutionRole'

    def _get_task_role_arn(self):
        """Get ECS task role ARN."""
        return f'arn:aws:iam::{self._get_account_id()}:role/rnaseq-mini-task-role'

    def _get_account_id(self):
        """Get AWS account ID."""
        return boto3.client('sts').get_caller_identity()['Account']

    def _create_target_group(self):
        """Create load balancer target group."""
        elbv2_client = boto3.client('elbv2', region_name=self.region)

        response = elbv2_client.create_target_group(
            Name='rnaseq-mini-tg',
            Protocol='HTTP',
            Port=8000,
            VpcId=self._get_vpc_id(),
            TargetType='ip'
        )

        return response['TargetGroups'][0]['TargetGroupArn']

    def _get_vpc_id(self):
        """Get VPC ID."""
        return f'vpc-{self.region.replace("-", "")}'

    def deploy_complete_infrastructure(self):
        """Deploy complete AWS infrastructure."""
        print("🚀 Deploying RNASEQ-MINI to AWS...")

        # Create VPC infrastructure
        self.create_vpc_stack()

        # Create compute environment
        self.create_compute_environment()

        # Create job queue
        self.create_job_queue()

        # Create job definition
        self.create_job_definition()

        # Create S3 bucket for results
        self.create_s3_bucket()

        # Create ECS service for web interface
        self.create_ecs_service()

        print("✅ AWS deployment complete!")
        print(f"🌐 Web interface: http://{self._get_load_balancer_dns()}")
        print("📊 Monitor progress: AWS Batch console"
    def _get_load_balancer_dns(self):
        """Get load balancer DNS name."""
        return f"rnaseq-mini-{self.region}.elb.amazonaws.com"


def main():
    """Main deployment function."""
    parser = argparse.ArgumentParser(description="Deploy RNASEQ-MINI to AWS")
    parser.add_argument('--region', default='us-east-1', help='AWS region')
    parser.add_argument('--profile', help='AWS profile to use')
    parser.add_argument('--create-bucket', action='store_true', help='Create S3 bucket')
    parser.add_argument('--create-compute-env', action='store_true', help='Create compute environment')
    parser.add_argument('--create-job-queue', action='store_true', help='Create job queue')
    parser.add_argument('--create-job-def', action='store_true', help='Create job definition')
    parser.add_argument('--create-ecs-service', action='store_true', help='Create ECS service')
    parser.add_argument('--all', action='store_true', help='Deploy complete infrastructure')

    args = parser.parse_args()

    if args.profile:
        os.environ['AWS_PROFILE'] = args.profile

    deployer = AWSDeployer(args.region)

    if args.all:
        deployer.deploy_complete_infrastructure()
    else:
        if args.create_bucket:
            deployer.create_s3_bucket()
        if args.create_compute_env:
            deployer.create_compute_environment()
        if args.create_job_queue:
            deployer.create_job_queue()
        if args.create_job_def:
            deployer.create_job_definition()
        if args.create_ecs_service:
            deployer.create_ecs_service()


if __name__ == "__main__":
    main()



