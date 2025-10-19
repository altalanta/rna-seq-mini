#!/usr/bin/env python3
"""
Snakemake wrapper for intelligent error handling.
Automatically classifies and recovers from errors in pipeline execution.
"""

import os
import sys
import json
import logging
import subprocess
import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any

# Import error handling components
from error_handler import (
    ErrorClassifier, ErrorContext, ErrorRecovery, ErrorReporter,
    ClassifiedError, ErrorSolution
)

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class SnakemakeErrorHandler:
    """Handles errors in Snakemake rule execution with intelligent recovery."""

    def __init__(self, rule_name: str, sample_id: str = None):
        self.rule_name = rule_name
        self.sample_id = sample_id or "unknown"
        self.classifier = ErrorClassifier()
        self.recovery = ErrorRecovery()
        self.reporter = ErrorReporter()

    def execute_with_recovery(self, command: str, max_retries: int = 3) -> bool:
        """Execute a command with intelligent error recovery."""
        logger.info(f"Executing rule '{self.rule_name}' for sample '{self.sample_id}'")

        for attempt in range(max_retries + 1):
            try:
                logger.info(f"Attempt {attempt + 1}/{max_retries + 1}: {command}")

                # Execute the command
                result = subprocess.run(
                    command,
                    shell=True,
                    capture_output=True,
                    text=True,
                    cwd=os.getcwd(),
                    timeout=3600  # 1 hour timeout
                )

                if result.returncode == 0:
                    logger.info(f"Rule '{self.rule_name}' completed successfully")
                    return True

                # Command failed, analyze the error
                error_message = result.stderr or result.stdout or "Unknown error"
                context = ErrorContext(
                    rule_name=self.rule_name,
                    sample_id=self.sample_id,
                    command=command,
                    exit_code=result.returncode,
                    stdout=result.stdout,
                    stderr=result.stderr,
                    working_directory=os.getcwd()
                )

                # Classify the error
                classified_error = self.classifier.classify_error(error_message, context)
                classified_error.retry_count = attempt

                # Record the error
                self.recovery.record_error(classified_error)

                # Check if we can retry
                if attempt < max_retries and classified_error.can_retry:
                    logger.warning(f"Rule '{self.rule_name}' failed (attempt {attempt + 1}), will retry...")
                    logger.warning(f"Error: {classified_error.root_cause}")

                    # Wait before retry (exponential backoff)
                    wait_time = min(30, 2 ** attempt)
                    logger.info(f"Waiting {wait_time} seconds before retry...")
                    import time
                    time.sleep(wait_time)

                    continue

                # Cannot retry or max retries reached
                logger.error(f"Rule '{self.rule_name}' failed after {attempt + 1} attempts")
                logger.error(f"Root cause: {classified_error.root_cause}")

                # Generate error report
                report_file = self.reporter.generate_report(classified_error)
                logger.error(f"Detailed error report: {report_file}")

                # Try to execute solutions
                self._try_automatic_solutions(classified_error)

                return False

            except subprocess.TimeoutExpired:
                logger.error(f"Rule '{self.rule_name}' timed out")
                if attempt < max_retries:
                    logger.info("Will retry with longer timeout...")
                    continue
                return False

            except KeyboardInterrupt:
                logger.info("Execution interrupted by user")
                return False

            except Exception as e:
                logger.error(f"Unexpected error in rule '{self.rule_name}': {e}")
                return False

    def _try_automatic_solutions(self, classified_error: ClassifiedError) -> bool:
        """Try to automatically execute error solutions."""
        logger.info("Attempting automatic error recovery...")

        for solution in classified_error.solutions:
            if solution.priority <= 2:  # Only try high-priority solutions automatically
                logger.info(f"Trying solution: {solution.title}")

                success = self.recovery.execute_solution(solution, classified_error.context)
                if success:
                    logger.info(f"Solution '{solution.title}' succeeded!")

                    # Mark error as resolved
                    self.recovery.mark_error_resolved(classified_error.error_id)

                    # Try the original command again
                    if self.execute_with_recovery(classified_error.context.command, 1):
                        logger.info("Error recovery successful!")
                        return True

        logger.warning("Automatic error recovery failed")
        return False

def main():
    """Main function for Snakemake integration."""
    import argparse

    parser = argparse.ArgumentParser(description='Snakemake error handling wrapper')
    parser.add_argument('--rule', required=True, help='Snakemake rule name')
    parser.add_argument('--sample', help='Sample ID')
    parser.add_argument('--command', required=True, help='Command to execute')
    parser.add_argument('--max-retries', type=int, default=3, help='Maximum retry attempts')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose logging')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Initialize error handler
    handler = SnakemakeErrorHandler(args.rule, args.sample)

    # Execute with error recovery
    success = handler.execute_with_recovery(args.command, args.max_retries)

    # Exit with appropriate code for Snakemake
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
