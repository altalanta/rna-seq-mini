#!/usr/bin/env python3
"""
Cost monitoring and optimization for RNASEQ-MINI cloud deployments.
Tracks AWS costs and provides recommendations for cost optimization.
"""

import os
import sys
import json
import logging
import argparse
import datetime
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
import boto3
from botocore.exceptions import ClientError, NoCredentialsError
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class CostMetrics:
    """Cost metrics for analysis."""
    service: str
    cost: float
    usage_quantity: float
    usage_unit: str
    resource_id: str
    tags: Dict[str, str]

@dataclass
class CostOptimization:
    """Cost optimization recommendations."""
    service: str
    current_cost: float
    potential_savings: float
    recommendation: str
    priority: str  # high, medium, low
    implementation_effort: str  # easy, medium, hard

class CostMonitor:
    """Monitors and analyzes AWS costs for RNASEQ-MINI deployments."""

    def __init__(self, region_name: str = 'us-east-1'):
        self.region_name = region_name

        try:
            self.cost_client = boto3.client('costexplorer', region_name=region_name)
            self.ce_client = boto3.client('cost-optimization-hub', region_name=region_name)
        except NoCredentialsError:
            logger.error("AWS credentials not found. Please configure AWS CLI or set environment variables.")
            sys.exit(1)

    def get_cost_data(self, start_date: str, end_date: str, granularity: str = 'DAILY') -> List[CostMetrics]:
        """Get cost data from AWS Cost Explorer."""
        try:
            # Get cost and usage data
            response = self.cost_client.get_cost_and_usage(
                TimePeriod={
                    'Start': start_date,
                    'End': end_date
                },
                Granularity=granularity,
                Metrics=['BlendedCost', 'UsageQuantity'],
                GroupBy=[
                    {
                        'Type': 'DIMENSION',
                        'Key': 'SERVICE'
                    },
                    {
                        'Type': 'DIMENSION',
                        'Key': 'RESOURCE_ID'
                    }
                ]
            )

            cost_data = []
            for result in response['ResultsByTime']:
                for group in result['Groups']:
                    service = None
                    resource_id = None

                    # Extract service and resource info
                    for key in group['Keys']:
                        if key.startswith('Service/'):
                            service = key.split('/')[-1]
                        elif key.startswith('ResourceId/'):
                            resource_id = key.split('/')[-1]

                    # Parse metrics
                    blended_cost = group['Metrics']['BlendedCost']['Amount']
                    usage_quantity = group['Metrics']['UsageQuantity']['Amount']
                    usage_unit = group['Metrics']['UsageQuantity']['Unit']

                    cost_data.append(CostMetrics(
                        service=service or 'Unknown',
                        cost=float(blended_cost),
                        usage_quantity=float(usage_quantity),
                        usage_unit=usage_unit,
                        resource_id=resource_id or 'N/A',
                        tags={}
                    ))

            return cost_data

        except ClientError as e:
            logger.error(f"Error getting cost data: {e}")
            return []

    def analyze_cost_patterns(self, cost_data: List[CostMetrics]) -> Dict[str, Any]:
        """Analyze cost patterns and identify trends."""
        if not cost_data:
            return {}

        # Convert to DataFrame for analysis
        df = pd.DataFrame([
            {
                'service': metric.service,
                'cost': metric.cost,
                'usage_quantity': metric.usage_quantity,
                'usage_unit': metric.usage_unit,
                'resource_id': metric.resource_id
            }
            for metric in cost_data
        ])

        # Group by service
        service_costs = df.groupby('service')['cost'].agg(['sum', 'mean', 'std', 'count']).round(2)

        # Calculate total cost
        total_cost = df['cost'].sum()

        # Identify top cost drivers
        top_services = service_costs.nlargest(5, 'sum')

        # Calculate cost distribution
        cost_distribution = (service_costs['sum'] / total_cost * 100).round(2)

        return {
            'total_cost': total_cost,
            'service_breakdown': service_costs.to_dict('index'),
            'top_cost_drivers': top_services.to_dict('index'),
            'cost_distribution': cost_distribution.to_dict(),
            'num_resources': len(df['resource_id'].unique())
        }

    def generate_optimization_recommendations(self, cost_data: List[CostMetrics]) -> List[CostOptimization]:
        """Generate cost optimization recommendations."""
        recommendations = []

        if not cost_data:
            return recommendations

        # Convert to DataFrame
        df = pd.DataFrame([
            {
                'service': metric.service,
                'cost': metric.cost,
                'usage_quantity': metric.usage_quantity,
                'usage_unit': metric.usage_unit
            }
            for metric in cost_data
        ])

        # Analyze each service for optimization opportunities
        service_totals = df.groupby('service')['cost'].sum()

        for service, total_cost in service_totals.items():
            service_data = df[df['service'] == service]

            # EC2 optimization
            if service == 'Amazon Elastic Compute Cloud - Compute':
                # Check for underutilized instances
                if total_cost > 100:  # Only recommend for significant costs
                    recommendations.append(CostOptimization(
                        service=service,
                        current_cost=total_cost,
                        potential_savings=total_cost * 0.3,  # Assume 30% potential savings
                        recommendation="Consider using Reserved Instances or Savings Plans for predictable workloads",
                        priority="high",
                        implementation_effort="medium"
                    ))

            # Batch optimization
            elif service == 'AWS Batch':
                if total_cost > 50:
                    recommendations.append(CostOptimization(
                        service=service,
                        current_cost=total_cost,
                        potential_savings=total_cost * 0.25,
                        recommendation="Implement auto-scaling to reduce idle compute time",
                        priority="high",
                        implementation_effort="easy"
                    ))

            # Storage optimization
            elif 'Storage' in service or 'S3' in service:
                if total_cost > 20:
                    recommendations.append(CostOptimization(
                        service=service,
                        current_cost=total_cost,
                        potential_savings=total_cost * 0.2,
                        recommendation="Use lifecycle policies and intelligent tiering for storage optimization",
                        priority="medium",
                        implementation_effort="easy"
                    ))

        # General recommendations
        if service_totals.sum() > 500:
            recommendations.append(CostOptimization(
                service="All Services",
                current_cost=service_totals.sum(),
                potential_savings=service_totals.sum() * 0.15,
                recommendation="Set up cost budgets and alerts in AWS Budgets",
                priority="medium",
                implementation_effort="easy"
            ))

        return recommendations

    def create_cost_report(self, cost_data: List[CostMetrics], output_file: str = None) -> Dict[str, Any]:
        """Create a comprehensive cost report."""
        # Analyze cost patterns
        analysis = self.analyze_cost_patterns(cost_data)

        # Generate recommendations
        recommendations = self.generate_optimization_recommendations(cost_data)

        # Calculate potential savings
        total_potential_savings = sum(rec.potential_savings for rec in recommendations)

        # Create report
        report = {
            'report_date': datetime.datetime.now().isoformat(),
            'analysis_period': 'Last 30 days',  # This would be dynamic
            'total_cost': analysis.get('total_cost', 0),
            'potential_savings': total_potential_savings,
            'savings_percentage': (total_potential_savings / max(analysis.get('total_cost', 1), 1)) * 100,
            'cost_analysis': analysis,
            'recommendations': [
                {
                    'service': rec.service,
                    'current_cost': rec.current_cost,
                    'potential_savings': rec.potential_savings,
                    'recommendation': rec.recommendation,
                    'priority': rec.priority,
                    'implementation_effort': rec.implementation_effort
                }
                for rec in recommendations
            ]
        }

        # Save report if output file specified
        if output_file:
            with open(output_file, 'w') as f:
                json.dump(report, f, indent=2)
            logger.info(f"Cost report saved to {output_file}")

        return report

    def plot_cost_trends(self, cost_data: List[CostMetrics], output_file: str = None):
        """Create cost trend visualizations."""
        if not cost_data:
            logger.warning("No cost data available for plotting")
            return

        # Convert to DataFrame
        df = pd.DataFrame([
            {
                'service': metric.service,
                'cost': metric.cost,
                'usage_quantity': metric.usage_quantity,
                'date': datetime.datetime.now()  # This would be actual dates from API
            }
            for metric in cost_data
        ])

        # Create plots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))

        # Cost by service (pie chart)
        service_costs = df.groupby('service')['cost'].sum().sort_values(ascending=False)
        axes[0, 0].pie(service_costs.values[:10], labels=service_costs.index[:10], autopct='%1.1f%%')
        axes[0, 0].set_title('Cost Distribution by Service')

        # Cost by service (bar chart)
        service_costs.head(10).plot(kind='bar', ax=axes[0, 1])
        axes[0, 1].set_title('Top 10 Cost Drivers')
        axes[0, 1].set_ylabel('Cost ($)')

        # Usage vs Cost correlation
        usage_by_service = df.groupby('service')['usage_quantity'].sum()
        cost_by_service = df.groupby('service')['cost'].sum()

        # Create scatter plot for top services
        top_services = cost_by_service.nlargest(10).index
        scatter_data = pd.DataFrame({
            'usage': usage_by_service[top_services],
            'cost': cost_by_service[top_services]
        })

        axes[1, 0].scatter(scatter_data['usage'], scatter_data['cost'])
        for i, service in enumerate(top_services):
            axes[1, 0].annotate(service, (scatter_data['usage'].iloc[i], scatter_data['cost'].iloc[i]))
        axes[1, 0].set_xlabel('Usage Quantity')
        axes[1, 0].set_ylabel('Cost ($)')
        axes[1, 0].set_title('Usage vs Cost Correlation')

        # Cost efficiency (cost per unit usage)
        efficiency = (cost_by_service / usage_by_service.replace(0, 1)).sort_values(ascending=False)
        efficiency.head(10).plot(kind='bar', ax=axes[1, 1])
        axes[1, 1].set_title('Cost Efficiency (Cost per Usage Unit)')
        axes[1, 1].set_ylabel('Cost per Unit')

        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Cost plots saved to {output_file}")
        else:
            plt.show()

def setup_cost_monitoring_dashboard():
    """Setup CloudWatch dashboards for cost monitoring."""
    try:
        # This would create CloudWatch dashboards for cost monitoring
        # For now, just log what would be created
        logger.info("Setting up cost monitoring dashboard...")

        dashboard_config = {
            'dashboardName': 'RNASEQ-MINI-Cost-Monitoring',
            'widgets': [
                {
                    'type': 'cost',
                    'properties': {
                        'metrics': ['BlendedCost'],
                        'groupBy': ['ServiceName'],
                        'period': 86400
                    }
                }
            ]
        }

        logger.info("Dashboard configuration ready")
        logger.info("Dashboard would include:")
        logger.info("- Daily cost trends")
        logger.info("- Service cost breakdown")
        logger.info("- Budget vs actual spending")
        logger.info("- Cost optimization recommendations")

        return True

    except Exception as e:
        logger.error(f"Error setting up dashboard: {e}")
        return False

def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(description='Cost monitoring for RNASEQ-MINI')
    parser.add_argument('--start-date', default=None, help='Start date (YYYY-MM-DD)')
    parser.add_argument('--end-date', default=None, help='End date (YYYY-MM-DD)')
    parser.add_argument('--region', default='us-east-1', help='AWS region')
    parser.add_argument('--output-report', help='Output file for cost report')
    parser.add_argument('--output-plot', help='Output file for cost plots')
    parser.add_argument('--setup-dashboard', action='store_true', help='Setup cost monitoring dashboard')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose logging')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Setup dashboard if requested
    if args.setup_dashboard:
        success = setup_cost_monitoring_dashboard()
        if not success:
            sys.exit(1)
        return

    # Default date range (last 30 days)
    if not args.start_date:
        args.start_date = (datetime.datetime.now() - datetime.timedelta(days=30)).strftime('%Y-%m-%d')
    if not args.end_date:
        args.end_date = datetime.datetime.now().strftime('%Y-%m-%d')

    logger.info(f"Analyzing costs from {args.start_date} to {args.end_date}")

    # Initialize cost monitor
    monitor = CostMonitor(region_name=args.region)

    # Get cost data
    cost_data = monitor.get_cost_data(args.start_date, args.end_date)

    if not cost_data:
        logger.warning("No cost data found for the specified period")
        return

    logger.info(f"Retrieved {len(cost_data)} cost records")

    # Create cost report
    report = monitor.create_cost_report(cost_data, args.output_report)

    # Print summary
    print("
💰 Cost Analysis Summary:"    print(f"   Total cost: ${report['total_cost']:.2f}")
    print(f"   Potential savings: ${report['potential_savings']:.2f} ({report['savings_percentage']:.1f}%)")
    print(f"   Number of services: {len(report['cost_analysis'].get('service_breakdown', {}))}")

    # Print top recommendations
    if report['recommendations']:
        print("
🎯 Top Cost Optimization Recommendations:"        for i, rec in enumerate(report['recommendations'][:5], 1):
            print(f"   {i}. {rec['service']}: Save ${rec['potential_savings']:.2f}")
            print(f"      {rec['recommendation']} ({rec['priority']} priority)")

    # Create plots
    if args.output_plot:
        monitor.plot_cost_trends(cost_data, args.output_plot)

    logger.info("Cost analysis completed")

if __name__ == "__main__":
    main()
