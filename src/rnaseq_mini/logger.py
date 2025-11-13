import logging
import logging.config
import os
from pathlib import Path

import structlog
import yaml

def setup_logging(config_path: Path = Path("config/logging.yaml"), log_format: str = "console"):
    """
    Sets up structured logging for the entire application.
    """
    if not config_path.exists():
        print(f"Warning: Logging config file not found at {config_path}. Using basic logging.")
        logging.basicConfig(level=logging.INFO)
        return

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Determine which formatter to use
    selected_formatter = log_format if log_format in config['formatters'] else 'console'
    
    # Configure standard logging from the dictionary
    logging.config.dictConfig(config)

    # Configure structlog to wrap the standard logger
    structlog.configure(
        processors=[
            structlog.stdlib.filter_by_level,
            structlog.stdlib.add_logger_name,
            structlog.stdlib.add_log_level,
            structlog.stdlib.PositionalArgumentsFormatter(),
            structlog.processors.TimeStamper(fmt="iso"),
            structlog.processors.StackInfoRenderer(),
            structlog.processors.format_exc_info,
            structlog.processors.UnicodeDecoder(),
            structlog.dev.ConsoleRenderer() if selected_formatter == 'console' else structlog.processors.JSONRenderer()
        ],
        context_class=dict,
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.stdlib.BoundLogger,
        cache_logger_on_first_use=True,
    )

def get_logger(name: str):
    """
    Returns a pre-configured logger instance.
    """
    return structlog.get_logger(name)

# --- Initial setup on import ---
# Use an environment variable to switch between formats
# e.g., `LOG_FORMAT=json python scripts/run_stage.py ...`
log_format_env = os.environ.get("LOG_FORMAT", "console")
setup_logging(log_format=log_format_env)




