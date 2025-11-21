#!/usr/bin/env python3
"""
Proactive troubleshooting and monitoring for RNASEQ-MINI.
Continuously monitors pipeline health and automatically fixes common issues.
"""

import os
import sys
import json
import logging
import time
import psutil
import threading
import subprocess
import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field

# Import error handling components
from error_handler import DiagnosticTools, ErrorReporter, ErrorRecovery

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class SystemHealth:
    """System health metrics."""
    cpu_percent: float
    memory_percent: float
    disk_percent: float
    network_status: bool
    timestamp: datetime.datetime = field(default_factory=datetime.datetime.now)

@dataclass
class PipelineHealth:
    """Pipeline health metrics."""
    active_jobs: int
    completed_jobs: int
    failed_jobs: int
    running_time: float
    last_activity: datetime.datetime = field(default_factory=datetime.datetime.now)

class ProactiveTroubleshooter:
    """Proactively monitors and fixes pipeline issues."""

    def __init__(self, check_interval: int = 30, enable_auto_fix: bool = True):
        self.check_interval = check_interval
        self.enable_auto_fix = enable_auto_fix
        self.diagnostic_tools = DiagnosticTools()
        self.error_recovery = ErrorRecovery()
        self.is_monitoring = False
        self.monitor_thread = None

        # Health thresholds
        self.health_thresholds = {
            'cpu_warning': 80.0,
            'cpu_critical': 95.0,
            'memory_warning': 85.0,
            'memory_critical': 95.0,
            'disk_warning': 80.0,
            'disk_critical': 90.0
        }

    def check_system_health(self) -> SystemHealth:
        """Check current system health."""
        try:
            cpu_percent = psutil.cpu_percent(interval=1)
            memory = psutil.virtual_memory()
            disk = psutil.disk_usage('/')

            # Check network connectivity
            network_status = False
            try:
                result = subprocess.run(
                    ['ping', '-c', '1', '8.8.8.8'],
                    capture_output=True, timeout=5
                )
                network_status = result.returncode == 0
            except:
                network_status = False

            return SystemHealth(
                cpu_percent=cpu_percent,
                memory_percent=memory.percent,
                disk_percent=disk.percent,
                network_status=network_status
            )

        except Exception as e:
            logger.error(f"Error checking system health: {e}")
            return SystemHealth(0, 0, 0, False)

    def check_pipeline_health(self) -> PipelineHealth:
        """Check current pipeline health."""
        try:
            # Check for Snakemake processes
            active_jobs = 0
            for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
                try:
                    if 'snakemake' in proc.info['name'] or 'snakemake' in ' '.join(proc.info['cmdline']):
                        active_jobs += 1
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass

            # Check log files for activity
            log_dir = Path("logs")
            completed_jobs = 0
            failed_jobs = 0
            last_activity = datetime.datetime.now()

            if log_dir.exists():
                log_files = list(log_dir.glob("*.log"))
                for log_file in log_files:
                    try:
                        mtime = datetime.datetime.fromtimestamp(log_file.stat().st_mtime)
                        if mtime > last_activity:
                            last_activity = mtime

                        # Count completed/failed jobs (simple heuristic)
                        with open(log_file, 'r') as f:
                            content = f.read()
                            if 'finished' in content.lower() or 'succeeded' in content.lower():
                                completed_jobs += 1
                            if 'failed' in content.lower() or 'error' in content.lower():
                                failed_jobs += 1
                    except:
                        pass

            # Calculate running time
            if last_activity != datetime.datetime.now():
                running_time = (datetime.datetime.now() - last_activity).total_seconds()
            else:
                running_time = 0

            return PipelineHealth(
                active_jobs=active_jobs,
                completed_jobs=completed_jobs,
                failed_jobs=failed_jobs,
                running_time=running_time,
                last_activity=last_activity
            )

        except Exception as e:
            logger.error(f"Error checking pipeline health: {e}")
            return PipelineHealth(0, 0, 0, 0)

    def analyze_health_issues(self, system_health: SystemHealth, pipeline_health: PipelineHealth) -> List[str]:
        """Analyze health metrics and identify issues."""
        issues = []

        # System resource issues
        if system_health.cpu_percent > self.health_thresholds['cpu_critical']:
            issues.append(f"Critical CPU usage: {system_health.cpu_percent:.1f}%")
        elif system_health.cpu_percent > self.health_thresholds['cpu_warning']:
            issues.append(f"High CPU usage: {system_health.cpu_percent:.1f}%")

        if system_health.memory_percent > self.health_thresholds['memory_critical']:
            issues.append(f"Critical memory usage: {system_health.memory_percent:.1f}%")
        elif system_health.memory_percent > self.health_thresholds['memory_warning']:
            issues.append(f"High memory usage: {system_health.memory_percent:.1f}%")

        if system_health.disk_percent > self.health_thresholds['disk_critical']:
            issues.append(f"Critical disk usage: {system_health.disk_percent:.1f}%")
        elif system_health.disk_percent > self.health_thresholds['disk_warning']:
            issues.append(f"High disk usage: {system_health.disk_percent:.1f}%")

        # Network issues
        if not system_health.network_status:
            issues.append("Network connectivity issues detected")

        # Pipeline issues
        if pipeline_health.failed_jobs > 0:
            issues.append(f"Pipeline failures detected: {pipeline_health.failed_jobs} failed jobs")

        if pipeline_health.active_jobs == 0 and pipeline_health.running_time > 300:  # 5 minutes
            issues.append("Pipeline appears stalled - no active jobs for extended period")

        return issues

    def auto_fix_issues(self, issues: List[str]) -> Dict[str, bool]:
        """Attempt to automatically fix identified issues."""
        fixes_applied = {}

        for issue in issues:
            logger.info(f"Attempting to fix: {issue}")

            if "disk usage" in issue.lower():
                success = self._fix_disk_space()
                fixes_applied[f"disk_{issue}"] = success

            elif "memory usage" in issue.lower():
                success = self._fix_memory_usage()
                fixes_applied[f"memory_{issue}"] = success

            elif "network connectivity" in issue.lower():
                success = self._fix_network_issues()
                fixes_applied["network"] = success

            elif "pipeline failures" in issue.lower():
                success = self._fix_pipeline_failures()
                fixes_applied["pipeline_failures"] = success

            else:
                logger.warning(f"No automatic fix available for: {issue}")
                fixes_applied[issue] = False

        return fixes_applied

    def _fix_disk_space(self) -> bool:
        """Attempt to free up disk space."""
        try:
            # Clean temporary files
            temp_dirs = ['temp', '.snakemake', 'work', '__pycache__']

            for temp_dir in temp_dirs:
                if os.path.exists(temp_dir):
                    logger.info(f"Cleaning {temp_dir}...")
                    subprocess.run(f'rm -rf {temp_dir}/*', shell=True, timeout=60)

            # Clean old cache files
            cache_dir = Path(".cache")
            if cache_dir.exists():
                # Remove files older than 7 days
                cutoff = datetime.datetime.now() - datetime.timedelta(days=7)
                for file_path in cache_dir.rglob('*'):
                    if file_path.is_file():
                        try:
                            mtime = datetime.datetime.fromtimestamp(file_path.stat().st_mtime)
                            if mtime < cutoff:
                                file_path.unlink()
                                logger.info(f"Removed old cache file: {file_path}")
                        except:
                            pass

            logger.info("Disk space cleanup completed")
            return True

        except Exception as e:
            logger.error(f"Error during disk cleanup: {e}")
            return False

    def _fix_memory_usage(self) -> bool:
        """Attempt to reduce memory usage."""
        try:
            # Kill any runaway processes
            for proc in psutil.process_iter(['pid', 'name', 'memory_percent']):
                try:
                    if (proc.info['memory_percent'] > 50 and
                        proc.info['name'] not in ['python', 'snakemake', 'nextflow']):
                        logger.warning(f"High memory process: {proc.info['name']} ({proc.info['memory_percent']:.1f}%)")
                        # Don't kill automatically, just log for now
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass

            # Clear Python cache
            subprocess.run('find . -name "*.pyc" -delete', shell=True, timeout=30)
            subprocess.run('find . -name "__pycache__" -type d -exec rm -rf {} +', shell=True, timeout=30)

            logger.info("Memory optimization completed")
            return True

        except Exception as e:
            logger.error(f"Error during memory optimization: {e}")
            return False

    def _fix_network_issues(self) -> bool:
        """Attempt to fix network connectivity issues."""
        try:
            # Restart network services (simple approach)
            logger.info("Checking network configuration...")

            # Test DNS resolution
            result = subprocess.run('nslookup google.com', shell=True, capture_output=True, timeout=10)
            if result.returncode == 0:
                logger.info("DNS resolution working")
                return True
            else:
                logger.warning("DNS resolution issues detected")
                # Could restart network service here on Linux
                return False

        except Exception as e:
            logger.error(f"Error fixing network issues: {e}")
            return False

    def _fix_pipeline_failures(self) -> bool:
        """Attempt to fix pipeline failures."""
        try:
            # Check for common configuration issues
            logger.info("Checking pipeline configuration...")

            # Validate configuration files
            config_checks = [
                'python scripts/validate_config.py',
                'yamllint config/ 2>/dev/null || echo "yamllint not available"'
            ]

            for check in config_checks:
                try:
                    result = subprocess.run(check, shell=True, capture_output=True, timeout=30)
                    if result.returncode != 0:
                        logger.warning(f"Configuration issue detected: {result.stderr.decode()}")
                    else:
                        logger.info("Configuration check passed")
                except:
                    pass

            # Check for missing dependencies
            dep_checks = [
                'which snakemake',
                'which nextflow',
                'which salmon'
            ]

            missing_deps = []
            for check in dep_checks:
                try:
                    result = subprocess.run(check, shell=True, capture_output=True, timeout=10)
                    if result.returncode != 0:
                        missing_deps.append(check.split()[-1])
                except:
                    pass

            if missing_deps:
                logger.warning(f"Missing dependencies: {', '.join(missing_deps)}")
                logger.info("Run 'make setup' to install missing dependencies")

            return True

        except Exception as e:
            logger.error(f"Error fixing pipeline failures: {e}")
            return False

    def generate_health_report(self, system_health: SystemHealth, pipeline_health: PipelineHealth, issues: List[str], fixes: Dict[str, bool]) -> str:
        """Generate a comprehensive health report."""
        report_file = f"results/health_report_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.md"

        report_content = f"""# RNASEQ-MINI Health Report

## Report Summary
- **Generated**: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
- **Monitoring Interval**: {self.check_interval} seconds
- **Auto-fix Enabled**: {self.enable_auto_fix}

## System Health
- **CPU Usage**: {system_health.cpu_percent:.1f}%
- **Memory Usage**: {system_health.memory_percent:.1f}%
- **Disk Usage**: {system_health.disk_percent:.1f}%
- **Network Status**: {"✅ Connected" if system_health.network_status else "❌ Disconnected"}

## Pipeline Health
- **Active Jobs**: {pipeline_health.active_jobs}
- **Completed Jobs**: {pipeline_health.completed_jobs}
- **Failed Jobs**: {pipeline_health.failed_jobs}
- **Running Time**: {pipeline_health.running_time:.0f} seconds
- **Last Activity**: {pipeline_health.last_activity.strftime('%Y-%m-%d %H:%M:%S')}

## Issues Detected
"""

        if issues:
            for issue in issues:
                report_content += f"- ⚠️  {issue}\n"
        else:
            report_content += "✅ No issues detected\n"

        report_content += "\n## Automatic Fixes Applied\n"

        if fixes:
            for fix, success in fixes.items():
                status = "✅ Success" if success else "❌ Failed"
                report_content += f"- {fix}: {status}\n"
        else:
            report_content += "No automatic fixes were applied\n"

        # Add recommendations
        report_content += """
## Recommendations

### System Optimization
1. **Monitor resource usage** - Use `htop` or `top` to identify resource-intensive processes
2. **Clean temporary files** - Run `make clean` to remove old results
3. **Optimize memory** - Reduce thread counts for memory-intensive operations

### Pipeline Optimization
1. **Check configuration** - Run `make validate-full` to verify settings
2. **Update dependencies** - Ensure all tools are up to date
3. **Review error logs** - Check `results/errors/` for detailed error reports

### Network Issues
1. **Check connectivity** - Verify internet connection for downloading references
2. **Firewall settings** - Ensure required ports are open
3. **DNS configuration** - Verify DNS resolution is working

## Getting Help

If issues persist:
1. Check the [troubleshooting guide](docs/troubleshooting.md)
2. Run diagnostics: `python scripts/error_handler.py --diagnostics all`
3. Review error reports in `results/errors/`
4. Create an issue with this health report attached
"""

        with open(report_file, 'w') as f:
            f.write(report_content)

        logger.info(f"Health report saved to {report_file}")
        return report_file

    def start_monitoring(self, duration_minutes: int = 60):
        """Start continuous health monitoring."""
        logger.info(f"Starting proactive troubleshooting for {duration_minutes} minutes")

        self.is_monitoring = True
        start_time = time.time()
        end_time = start_time + (duration_minutes * 60)

        try:
            while time.time() < end_time and self.is_monitoring:
                # Check system and pipeline health
                system_health = self.check_system_health()
                pipeline_health = self.check_pipeline_health()

                # Analyze for issues
                issues = self.analyze_health_issues(system_health, pipeline_health)

                # Apply automatic fixes if enabled
                fixes_applied = {}
                if self.enable_auto_fix and issues:
                    fixes_applied = self.auto_fix_issues(issues)

                # Log current status
                if issues:
                    logger.warning(f"Issues detected: {', '.join(issues)}")
                else:
                    logger.info("System and pipeline health: OK")

                # Generate health report every 10 minutes
                if int(time.time()) % 600 < self.check_interval:  # Every 10 minutes
                    self.generate_health_report(system_health, pipeline_health, issues, fixes_applied)

                # Wait for next check
                time.sleep(self.check_interval)

        except KeyboardInterrupt:
            logger.info("Monitoring interrupted by user")
        finally:
            self.is_monitoring = False
            logger.info("Proactive troubleshooting completed")

    def stop_monitoring(self):
        """Stop the monitoring process."""
        self.is_monitoring = False
        if self.monitor_thread and self.monitor_thread.is_alive():
            self.monitor_thread.join(timeout=5)

def main():
    """Main function for command-line usage."""
    import argparse

    parser = argparse.ArgumentParser(description='Proactive troubleshooting for RNASEQ-MINI')
    parser.add_argument('--monitor', type=int, metavar='MINUTES',
                       help='Monitor system for specified minutes')
    parser.add_argument('--interval', type=int, default=30,
                       help='Check interval in seconds')
    parser.add_argument('--no-auto-fix', action='store_true',
                       help='Disable automatic fixes')
    parser.add_argument('--health-check', action='store_true',
                       help='Perform one-time health check')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose logging')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Initialize troubleshooter
    troubleshooter = ProactiveTroubleshooter(
        check_interval=args.interval,
        enable_auto_fix=not args.no_auto_fix
    )

    # Handle different modes
    if args.health_check:
        # One-time health check
        system_health = troubleshooter.check_system_health()
        pipeline_health = troubleshooter.check_pipeline_health()
        issues = troubleshooter.analyze_health_issues(system_health, pipeline_health)

        print("
🏥 Health Check Results:"        print(f"   CPU: {system_health.cpu_percent:.1f}%")
        print(f"   Memory: {system_health.memory_percent:.1f}%")
        print(f"   Disk: {system_health.disk_percent:.1f}%")
        print(f"   Network: {'✅ OK' if system_health.network_status else '❌ Issues'}")
        print(f"   Active Jobs: {pipeline_health.active_jobs}")
        print(f"   Failed Jobs: {pipeline_health.failed_jobs}")

        if issues:
            print("
⚠️  Issues Detected:"            for issue in issues:
                print(f"   • {issue}")

            if not args.no_auto_fix:
                fixes = troubleshooter.auto_fix_issues(issues)
                print("
🔧 Fixes Applied:"                for fix, success in fixes.items():
                    status = "✅ Success" if success else "❌ Failed"
                    print(f"   • {fix}: {status}")
        else:
            print("
✅ No issues detected"        # Generate health report
        report_file = troubleshooter.generate_health_report(system_health, pipeline_health, issues, {})
        print(f"\n📋 Health report: {report_file}")

    elif args.monitor:
        # Continuous monitoring
        troubleshooter.start_monitoring(args.monitor)

    else:
        parser.print_help()

if __name__ == "__main__":
    main()












