"""
Logging configuration for NeoSWGA.

Provides consistent logging across all modules with configurable levels and formats.
"""

import logging
import sys
from pathlib import Path
from typing import Optional
import os


# Default log format
DEFAULT_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
DETAILED_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s"
SIMPLE_FORMAT = "%(levelname)s: %(message)s"

# Color codes for terminal output
COLORS = {
    'DEBUG': '\033[36m',     # Cyan
    'INFO': '\033[32m',      # Green
    'WARNING': '\033[33m',   # Yellow
    'ERROR': '\033[31m',     # Red
    'CRITICAL': '\033[35m',  # Magenta
    'RESET': '\033[0m'       # Reset
}


class ColoredFormatter(logging.Formatter):
    """Custom formatter with color support for terminal output."""

    def format(self, record):
        if sys.stdout.isatty():
            # Add color to levelname
            levelname = record.levelname
            if levelname in COLORS:
                record.levelname = f"{COLORS[levelname]}{levelname}{COLORS['RESET']}"

        return super().format(record)


def setup_logging(
    level: str = "INFO",
    log_file: Optional[str] = None,
    format_style: str = "default",
    module_levels: Optional[dict] = None
):
    """
    Setup logging configuration for NeoSWGA.

    Args:
        level: Global logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional file path to write logs to
        format_style: Format style ('default', 'detailed', 'simple')
        module_levels: Optional dict of module-specific levels

    Example:
        >>> from neoswga.logging_config import setup_logging
        >>> setup_logging(level='DEBUG', log_file='neoswga.log')
        >>> setup_logging(level='INFO', module_levels={'neoswga.thermodynamics': 'DEBUG'})
    """

    # Choose format
    if format_style == "detailed":
        log_format = DETAILED_FORMAT
    elif format_style == "simple":
        log_format = SIMPLE_FORMAT
    else:
        log_format = DEFAULT_FORMAT

    # Get root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, level.upper()))

    # Remove existing handlers
    root_logger.handlers = []

    # Console handler with colors
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(getattr(logging, level.upper()))
    console_formatter = ColoredFormatter(log_format)
    console_handler.setFormatter(console_formatter)
    root_logger.addHandler(console_handler)

    # File handler (if specified)
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)  # Always log everything to file
        file_formatter = logging.Formatter(DETAILED_FORMAT)
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)

    # Set module-specific levels
    if module_levels:
        for module_name, module_level in module_levels.items():
            module_logger = logging.getLogger(module_name)
            module_logger.setLevel(getattr(logging, module_level.upper()))


def get_logger(name: str) -> logging.Logger:
    """
    Get logger for a module.

    Args:
        name: Module name (typically __name__)

    Returns:
        Configured logger

    Example:
        >>> from neoswga.logging_config import get_logger
        >>> logger = get_logger(__name__)
        >>> logger.info("Starting calculation...")
    """
    return logging.getLogger(name)


def configure_from_env():
    """
    Configure logging from environment variables.

    Environment variables:
        NeoSWGA_LOG_LEVEL: Logging level (default: INFO)
        NeoSWGA_LOG_FILE: Log file path (optional)
        NeoSWGA_LOG_FORMAT: Format style (default, detailed, simple)

    Example:
        $ export NeoSWGA_LOG_LEVEL=DEBUG
        $ export NeoSWGA_LOG_FILE=neoswga.log
        $ python -m neoswga ...
    """
    level = os.environ.get('NeoSWGA_LOG_LEVEL', 'INFO')
    log_file = os.environ.get('NeoSWGA_LOG_FILE', None)
    format_style = os.environ.get('NeoSWGA_LOG_FORMAT', 'default')

    setup_logging(level=level, log_file=log_file, format_style=format_style)


# Context manager for temporary log level
class temporary_log_level:
    """
    Context manager to temporarily change log level.

    Example:
        >>> from neoswga.logging_config import get_logger, temporary_log_level
        >>> logger = get_logger(__name__)
        >>> with temporary_log_level('DEBUG'):
        ...     logger.debug("This will be shown")
        >>> logger.debug("This won't be shown if level was INFO")
    """

    def __init__(self, level: str):
        self.level = level
        self.old_level = None

    def __enter__(self):
        logger = logging.getLogger()
        self.old_level = logger.level
        logger.setLevel(getattr(logging, self.level.upper()))
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        logger = logging.getLogger()
        logger.setLevel(self.old_level)


# Progress logging helper
class ProgressLogger:
    """
    Helper for logging progress of long-running operations.

    Example:
        >>> from neoswga.logging_config import ProgressLogger, get_logger
        >>> logger = get_logger(__name__)
        >>> progress = ProgressLogger(logger, total=1000, step=100)
        >>> for i in range(1000):
        ...     # Do work
        ...     progress.update(1)
    """

    def __init__(self, logger: logging.Logger, total: int, step: int = 10,
                 message: str = "Progress"):
        """
        Initialize progress logger.

        Args:
            logger: Logger to use
            total: Total number of items
            step: Log every N% progress
            message: Progress message
        """
        self.logger = logger
        self.total = total
        self.step = step
        self.message = message
        self.current = 0
        self.last_logged_pct = 0

    def update(self, n: int = 1):
        """Update progress by n items."""
        self.current += n
        pct = int((self.current / self.total) * 100)

        # Log at step intervals
        if pct >= self.last_logged_pct + self.step:
            self.logger.info(f"{self.message}: {pct}% ({self.current}/{self.total})")
            self.last_logged_pct = pct

    def finish(self):
        """Log completion."""
        self.logger.info(f"{self.message}: Complete ({self.total}/{self.total})")


# Initialize default logging
def init_default_logging():
    """Initialize default logging configuration."""
    # Check if already initialized
    if logging.getLogger().handlers:
        return

    # Try to load from environment
    if any(k.startswith('NeoSWGA_LOG_') for k in os.environ):
        configure_from_env()
    else:
        # Use defaults
        setup_logging(level='INFO', format_style='default')


# Auto-initialize on import
init_default_logging()


if __name__ == "__main__":
    # Demo logging configuration
    print("NeoSWGA Logging Configuration Demo")
    print("=" * 60)

    # Default setup
    print("\n1. Default logging (INFO level):")
    setup_logging(level='INFO')
    logger = get_logger(__name__)
    logger.debug("This won't be shown")
    logger.info("This will be shown")
    logger.warning("This is a warning")
    logger.error("This is an error")

    # Debug level
    print("\n2. Debug logging:")
    setup_logging(level='DEBUG')
    logger.debug("Now debug messages are shown")

    # With file output
    print("\n3. Logging to file:")
    setup_logging(level='INFO', log_file='/tmp/neoswga_test.log')
    logger.info("This goes to both console and file")
    print("Check /tmp/neoswga_test.log")

    # Module-specific levels
    print("\n4. Module-specific levels:")
    setup_logging(
        level='INFO',
        module_levels={
            'neoswga.thermodynamics': 'DEBUG',
            'neoswga.genetic_algorithm': 'WARNING'
        }
    )
    logger.info("Root logger at INFO")
    thermo_logger = get_logger('neoswga.thermodynamics')
    thermo_logger.debug("Thermodynamics at DEBUG")

    # Progress logging
    print("\n5. Progress logging:")
    progress = ProgressLogger(logger, total=100, step=25)
    for i in range(100):
        progress.update(1)
    progress.finish()

    # Temporary log level
    print("\n6. Temporary log level:")
    setup_logging(level='INFO')
    logger.debug("This won't be shown")
    with temporary_log_level('DEBUG'):
        logger.debug("This will be shown temporarily")
    logger.debug("This won't be shown again")

    print("\n" + "=" * 60)
    print("Demo complete!")
