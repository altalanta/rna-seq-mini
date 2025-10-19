#!/usr/bin/env python3
"""
Intelligent error handling and automated troubleshooting for RNASEQ-MINI.
Provides comprehensive error classification, recovery, and user-friendly reporting.
"""

import os
import sys
import json
import logging
import traceback
import subprocess
import datetime
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Callable
from dataclasses import dataclass, field
from enum import Enum
import time

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ErrorSeverity(Enum):
    """Error severity levels."""
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"
    CRITICAL = "critical"

class ErrorCategory(Enum):
    """Error categories for classification."""
    SYSTEM = "system"
    DEPENDENCY = "dependency"
    CONFIGURATION = "configuration"
    DATA = "data"
    COMPUTATIONAL = "computational"
    NETWORK = "network"
    PERMISSION = "permission"
    RESOURCE = "resource"
    UNKNOWN = "unknown"

@dataclass
class ErrorContext:
    """Context information for error analysis."""
    rule_name: str = ""
    sample_id: str = ""
    stage: str = ""
    command: str = ""
    exit_code: int = 0
    stdout: str = ""
    stderr: str = ""
    working_directory: str = ""
    environment: Dict[str, str] = field(default_factory=dict)
    timestamp: datetime.datetime = field(default_factory=datetime.datetime.now)

@dataclass
class ErrorSolution:
    """Solution information for errors."""
    title: str
    description: str
    commands: List[str]
    priority: int
    estimated_time: str
    success_probability: float

@dataclass
class ClassifiedError:
    """Classified error with metadata and solutions."""
    error_id: str
    original_message: str
    category: ErrorCategory
    severity: ErrorSeverity
    context: ErrorContext
    root_cause: str
    solutions: List[ErrorSolution]
    retry_count: int = 0
    max_retries: int = 3
    can_retry: bool = True

class ErrorClassifier:
    """Intelligent error classification system."""

    def __init__(self):
        # Error pattern definitions
        self.error_patterns = {
            ErrorCategory.SYSTEM: [
                (r"command not found", "Missing system command"),
                (r"no such file or directory", "File or directory not found"),
                (r"permission denied", "Insufficient permissions"),
                (r"disk.*full", "Insufficient disk space"),
                (r"out of memory", "Memory allocation failure"),
            ],
            ErrorCategory.DEPENDENCY: [
                (r"module.*not found", "Missing Python module"),
                (r"import.*failed", "Import error"),
                (r"conda.*environment", "Conda environment issue"),
                (r"pip.*install", "Package installation failure"),
            ],
            ErrorCategory.CONFIGURATION: [
                (r"invalid.*parameter", "Configuration parameter error"),
                (r"missing.*parameter", "Missing required parameter"),
                (r"yaml.*syntax", "YAML configuration syntax error"),
                (r"json.*syntax", "JSON configuration syntax error"),
            ],
            ErrorCategory.DATA: [
                (r"corrupt.*file", "Corrupted input file"),
                (r"format.*error", "File format error"),
                (r"empty.*file", "Empty input file"),
                (r"invalid.*fastq", "Invalid FASTQ format"),
            ],
            ErrorCategory.COMPUTATIONAL: [
                (r"numerical.*instability", "Numerical computation error"),
                (r"matrix.*singular", "Linear algebra error"),
                (r"convergence.*failed", "Algorithm convergence failure"),
                (r"optimization.*failed", "Optimization failure"),
            ],
            ErrorCategory.NETWORK: [
                (r"connection.*refused", "Network connection failure"),
                (r"timeout", "Network timeout"),
                (r"dns.*resolution", "DNS resolution failure"),
                (r"ssl.*certificate", "SSL/TLS certificate error"),
            ],
            ErrorCategory.RESOURCE: [
                (r"cpu.*limit", "CPU resource limit exceeded"),
                (r"memory.*limit", "Memory resource limit exceeded"),
                (r"disk.*quota", "Disk quota exceeded"),
                (r"job.*queue", "Job queue full"),
            ],
        }

        # Solution database
        self.solution_database = self._build_solution_database()

    def _build_solution_database(self) -> Dict[str, List[ErrorSolution]]:
        """Build database of error solutions."""
        return {
            "command_not_found": [
                ErrorSolution(
                    title="Install missing command",
                    description="The required command is not installed on the system",
                    commands=[
                        "conda install -c conda-forge {command_name}",
                        "apt-get install {command_name}  # On Ubuntu/Debian",
                        "yum install {command_name}  # On CentOS/RHEL"
                    ],
                    priority=1,
                    estimated_time="5 minutes",
                    success_probability=0.9
                )
            ],
            "permission_denied": [
                ErrorSolution(
                    title="Fix file permissions",
                    description="Insufficient permissions to access files or directories",
                    commands=[
                        "chmod +rwx {file_path}",
                        "chown {user}:{group} {file_path}",
                        "sudo chmod 755 {directory_path}"
                    ],
                    priority=1,
                    estimated_time="2 minutes",
                    success_probability=0.95
                )
            ],
            "disk_full": [
                ErrorSolution(
                    title="Free up disk space",
                    description="Insufficient disk space for analysis",
                    commands=[
                        "df -h  # Check disk usage",
                        "du -sh * | sort -hr | head -10  # Find large files",
                        "rm -rf results/temp/*  # Clean temporary files",
                        "find . -name "*.tmp" -delete  # Remove temp files"
                    ],
                    priority=1,
                    estimated_time="10 minutes",
                    success_probability=0.8
                )
            ],
            "out_of_memory": [
                ErrorSolution(
                    title="Reduce memory usage",
                    description="Process exceeded available memory",
                    commands=[
                        "Reduce thread count in config",
                        "Process samples in smaller batches",
                        "Use memory-efficient algorithms",
                        "Add swap space: sudo fallocate -l 8G /swapfile"
                    ],
                    priority=2,
                    estimated_time="15 minutes",
                    success_probability=0.7
                )
            ],
            "conda_environment": [
                ErrorSolution(
                    title="Fix conda environment",
                    description="Conda environment configuration issue",
                    commands=[
                        "conda env list  # Check environments",
                        "conda activate rnaseq-mini-base",
                        "mamba env create -f envs/base.yml  # Recreate environment",
                        "conda clean --all  # Clean conda cache"
                    ],
                    priority=2,
                    estimated_time="20 minutes",
                    success_probability=0.85
                )
            ],
            "configuration_error": [
                ErrorSolution(
                    title="Fix configuration",
                    description="Configuration file has errors",
                    commands=[
                        "python scripts/validate_config.py",
                        "Check YAML syntax with: yamllint config/params.yaml",
                        "Compare with example: diff config/params.yaml config/params_example.yaml"
                    ],
                    priority=1,
                    estimated_time="10 minutes",
                    success_probability=0.9
                )
            ]
        }

    def classify_error(self, error_message: str, context: ErrorContext) -> ClassifiedError:
        """Classify an error and determine solutions."""
        error_id = f"ERR_{int(time.time() * 1000)}"

        # Find matching patterns
        category = ErrorCategory.UNKNOWN
        severity = ErrorSeverity.MEDIUM
        root_cause = "Unknown error"
        solutions = []

        for cat, patterns in self.error_patterns.items():
            for pattern, description in patterns:
                if re.search(pattern, error_message, re.IGNORECASE):
                    category = cat
                    root_cause = description

                    # Determine severity based on category and message
                    if any(word in error_message.lower() for word in ['critical', 'fatal', 'segmentation fault']):
                        severity = ErrorSeverity.CRITICAL
                    elif cat in [ErrorCategory.SYSTEM, ErrorCategory.PERMISSION]:
                        severity = ErrorSeverity.HIGH
                    elif cat in [ErrorCategory.DATA, ErrorCategory.RESOURCE]:
                        severity = ErrorSeverity.MEDIUM

                    # Get solutions for this error type
                    solutions = self.solution_database.get(pattern.replace('\\', '').replace('.*', '').lower(), [])

                    break
            if category != ErrorCategory.UNKNOWN:
                break

        # Determine if error can be retried
        can_retry = category not in [ErrorCategory.CONFIGURATION, ErrorCategory.DATA]
        max_retries = 3 if can_retry else 0

        return ClassifiedError(
            error_id=error_id,
            original_message=error_message,
            category=category,
            severity=severity,
            context=context,
            root_cause=root_cause,
            solutions=solutions,
            max_retries=max_retries,
            can_retry=can_retry
        )

class ErrorRecovery:
    """Handles error recovery and retry logic."""

    def __init__(self, error_db_path: str = "results/.error_database.json"):
        self.error_db_path = Path(error_db_path)
        self.error_history = self._load_error_history()

    def _load_error_history(self) -> Dict[str, Any]:
        """Load error history from database."""
        if self.error_db_path.exists():
            try:
                with open(self.error_db_path, 'r') as f:
                    return json.load(f)
            except Exception as e:
                logger.warning(f"Could not load error history: {e}")
        return {"errors": [], "statistics": {}}

    def _save_error_history(self):
        """Save error history to database."""
        self.error_db_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.error_db_path, 'w') as f:
            json.dump(self.error_history, f, indent=2, default=str)

    def record_error(self, classified_error: ClassifiedError):
        """Record an error in the history database."""
        error_record = {
            "error_id": classified_error.error_id,
            "timestamp": classified_error.context.timestamp.isoformat(),
            "rule": classified_error.context.rule_name,
            "sample": classified_error.context.sample_id,
            "category": classified_error.category.value,
            "severity": classified_error.severity.value,
            "root_cause": classified_error.root_cause,
            "message": classified_error.original_message,
            "resolved": False,
            "retry_count": classified_error.retry_count
        }

        self.error_history["errors"].append(error_record)
        self._save_error_history()

    def can_retry_error(self, error_id: str) -> bool:
        """Check if an error can be retried."""
        for error_record in self.error_history["errors"]:
            if error_record["error_id"] == error_id:
                return error_record["retry_count"] < 3 and not error_record["resolved"]
        return False

    def mark_error_resolved(self, error_id: str):
        """Mark an error as resolved."""
        for error_record in self.error_history["errors"]:
            if error_record["error_id"] == error_id:
                error_record["resolved"] = True
                break
        self._save_error_history()

    def execute_solution(self, solution: ErrorSolution, context: ErrorContext) -> bool:
        """Execute a solution for an error."""
        logger.info(f"Executing solution: {solution.title}")

        for command_template in solution.commands:
            try:
                # Replace placeholders in command
                command = self._format_command(command_template, context)

                logger.info(f"Running: {command}")
                result = subprocess.run(command, shell=True, capture_output=True, text=True, timeout=300)

                if result.returncode == 0:
                    logger.info(f"Solution command succeeded: {command}")
                    return True
                else:
                    logger.warning(f"Solution command failed: {command}")
                    logger.warning(f"STDERR: {result.stderr}")

            except subprocess.TimeoutExpired:
                logger.error(f"Solution command timed out: {command}")
            except Exception as e:
                logger.error(f"Error executing solution: {e}")

        return False

    def _format_command(self, command_template: str, context: ErrorContext) -> str:
        """Format command template with context variables."""
        command = command_template

        # Extract file path from error message if possible
        file_match = re.search(r"'([^']+)'|\"([^\"]+)\"|(\S+\.\w+)", context.stderr + " " + context.stdout)
        if file_match:
            file_path = next((g for g in file_match.groups() if g), "")
            command = command.replace("{file_path}", file_path)

        # Extract directory path
        if "{directory_path}" in command and file_path:
            command = command.replace("{directory_path}", str(Path(file_path).parent))

        # Extract command name
        if "{command_name}" in command:
            cmd_match = re.search(r"(\w+): command not found", context.stderr)
            if cmd_match:
                command = command.replace("{command_name}", cmd_match.group(1))

        return command

class ErrorReporter:
    """Generates user-friendly error reports."""

    def __init__(self, output_dir: str = "results/errors"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_report(self, classified_error: ClassifiedError) -> str:
        """Generate a comprehensive error report."""
        timestamp = classified_error.context.timestamp.strftime("%Y%m%d_%H%M%S")
        report_file = self.output_dir / f"error_report_{classified_error.error_id}_{timestamp}.md"

        report_content = f"""# RNASEQ-MINI Error Report

## Error Summary
- **Error ID**: {classified_error.error_id}
- **Timestamp**: {classified_error.context.timestamp}
- **Category**: {classified_error.category.value}
- **Severity**: {classified_error.severity.value}
- **Rule**: {classified_error.context.rule_name or 'N/A'}
- **Sample**: {classified_error.context.sample_id or 'N/A'}

## Root Cause
{classified_error.root_cause}

## Original Error Message
```
{classified_error.original_message}
```

## Context Information
- **Working Directory**: {classified_error.context.working_directory}
- **Command**: {classified_error.context.command}
- **Exit Code**: {classified_error.context.exit_code}

## Standard Output
```
{classified_error.context.stdout[:1000]}{'...' if len(classified_error.context.stdout) > 1000 else ''}
```

## Standard Error
```
{classified_error.context.stderr[:1000]}{'...' if len(classified_error.context.stderr) > 1000 else ''}
```

## Recommended Solutions
"""

        for i, solution in enumerate(classified_error.solutions, 1):
            report_content += f"""
### Solution {i}: {solution.title}
**Priority**: {solution.priority} | **Time**: {solution.estimated_time} | **Success Rate**: {solution.success_probability:.1%}

**Description**: {solution.description}

**Commands to try**:
```bash
"""

            for cmd in solution.commands:
                report_content += f"{cmd}\n"

            report_content += "```\n\n"

        # Add troubleshooting tips
        report_content += """
## General Troubleshooting Tips

1. **Check system resources**: Ensure sufficient disk space, memory, and CPU cores
2. **Verify file permissions**: Make sure all input files are readable
3. **Validate configuration**: Run `make validate-full` to check configuration
4. **Check dependencies**: Verify all required tools are installed
5. **Review logs**: Check the full log files for additional details

## Getting Help

If these solutions don't resolve the issue:
1. Check the [troubleshooting guide](docs/troubleshooting.md)
2. Search existing [GitHub issues](https://github.com/rnaseq-mini/rnaseq-mini/issues)
3. Create a new issue with this error report attached
"""

        with open(report_file, 'w') as f:
            f.write(report_content)

        logger.info(f"Error report saved to {report_file}")
        return str(report_file)

class DiagnosticTools:
    """Automated diagnostic and troubleshooting tools."""

    def __init__(self):
        self.diagnostic_commands = {
            "system_check": [
                "Check system resources",
                ["df -h", "free -h", "nproc", "ulimit -a"]
            ],
            "dependency_check": [
                "Check software dependencies",
                ["conda env list", "python --version", "which snakemake", "which nextflow"]
            ],
            "configuration_check": [
                "Validate configuration files",
                ["python scripts/validate_config.py", "yamllint config/"]
            ],
            "data_check": [
                "Check input data integrity",
                ["ls -la config/samples.tsv", "head -5 config/samples.tsv"]
            ],
            "network_check": [
                "Check network connectivity",
                ["ping -c 3 google.com", "curl -I https://github.com"]
            ]
        }

    def run_diagnostics(self, category: str = "all") -> Dict[str, Any]:
        """Run diagnostic checks."""
        results = {}

        if category == "all":
            categories = self.diagnostic_commands.keys()
        else:
            categories = [category]

        for cat in categories:
            if cat in self.diagnostic_commands:
                title, commands = self.diagnostic_commands[cat]
                logger.info(f"Running {title}...")

                category_results = []
                for cmd in commands:
                    try:
                        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
                        category_results.append({
                            "command": cmd,
                            "success": result.returncode == 0,
                            "stdout": result.stdout.strip(),
                            "stderr": result.stderr.strip()
                        })
                    except subprocess.TimeoutExpired:
                        category_results.append({
                            "command": cmd,
                            "success": False,
                            "stdout": "",
                            "stderr": "Command timed out"
                        })
                    except Exception as e:
                        category_results.append({
                            "command": cmd,
                            "success": False,
                            "stdout": "",
                            "stderr": str(e)
                        })

                results[cat] = {
                    "title": title,
                    "commands": category_results,
                    "overall_success": any(r["success"] for r in category_results)
                }

        return results

def main():
    """Main function for command-line usage."""
    import argparse

    parser = argparse.ArgumentParser(description='Intelligent error handling for RNASEQ-MINI')
    parser.add_argument('--error-message', help='Error message to classify')
    parser.add_argument('--context-file', help='JSON file with error context')
    parser.add_argument('--diagnostics', choices=['system', 'dependency', 'configuration', 'data', 'network', 'all'],
                       default='all', help='Run diagnostic checks')
    parser.add_argument('--output-report', help='Output file for error report')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose logging')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Run diagnostics if requested
    if args.diagnostics:
        diagnostic_tools = DiagnosticTools()
        results = diagnostic_tools.run_diagnostics(args.diagnostics)

        print("
🔍 Diagnostic Results:"        for category, result in results.items():
            status = "✅ PASS" if result["overall_success"] else "❌ FAIL"
            print(f"   {category}: {status}")
            for cmd_result in result["commands"]:
                cmd_status = "✅" if cmd_result["success"] else "❌"
                print(f"     {cmd_status} {cmd_result['command']}")

        return

    # Classify error if provided
    if args.error_message:
        # Load context if provided
        context = ErrorContext()
        if args.context_file and os.path.exists(args.context_file):
            try:
                with open(args.context_file, 'r') as f:
                    context_data = json.load(f)
                    context = ErrorContext(**context_data)
            except Exception as e:
                logger.warning(f"Could not load context file: {e}")

        # Classify the error
        classifier = ErrorClassifier()
        classified_error = classifier.classify_error(args.error_message, context)

        # Generate report
        reporter = ErrorReporter()
        report_file = reporter.generate_report(classified_error)

        # Print summary
        print("
🚨 Error Analysis:"        print(f"   Category: {classified_error.category.value}")
        print(f"   Severity: {classified_error.severity.value}")
        print(f"   Root Cause: {classified_error.root_cause}")
        print(f"   Can Retry: {classified_error.can_retry}")
        print(f"   Solutions: {len(classified_error.solutions)}")

        if classified_error.solutions:
            print("
💡 Recommended Solutions:"            for i, solution in enumerate(classified_error.solutions[:3], 1):
                print(f"   {i}. {solution.title} (Priority: {solution.priority})")
                print(f"      {solution.description}")

        print(f"\n📋 Full report: {report_file}")

if __name__ == "__main__":
    main()
