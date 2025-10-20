"""
Structured logging infrastructure for Nommo Engine.

Provides consistent logging across all modules with support for:
- Structured JSON logging
- Performance metrics
- Simulation event tracking
- Debug tracing
"""

import json
import logging
import logging.handlers
import sys
import time
from datetime import datetime
from functools import wraps
from pathlib import Path
from typing import Any

from rich.console import Console
from rich.logging import RichHandler


class StructuredFormatter(logging.Formatter):
    """JSON formatter for structured logging."""

    def format(self, record: logging.LogRecord) -> str:
        log_data = {
            "timestamp": datetime.utcnow().isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
            "module": record.module,
            "function": record.funcName,
            "line": record.lineno,
        }

        if hasattr(record, "extra_data"):
            log_data.update(record.extra_data)

        if record.exc_info:
            log_data["exception"] = self.formatException(record.exc_info)

        return json.dumps(log_data)


class SimulationLogger:
    """Specialized logger for simulation events."""

    def __init__(self, name: str = "nommo.simulation"):
        self.logger = logging.getLogger(name)

    def log_event(
        self,
        event_type: str,
        universe_id: str,
        tick: int,
        data: dict[str, Any],
        level: int = logging.INFO,
    ) -> None:
        """Log a simulation event with structured data."""
        extra_data = {"event_type": event_type, "universe_id": universe_id, "tick": tick, **data}
        self.logger.log(level, f"Simulation event: {event_type}", extra={"extra_data": extra_data})

    def log_metrics(self, universe_id: str, tick: int, metrics: dict[str, Any]) -> None:
        """Log simulation metrics."""
        self.log_event("metrics", universe_id, tick, {"metrics": metrics})

    def log_emergence(
        self, universe_id: str, tick: int, emergence_type: str, details: dict[str, Any]
    ) -> None:
        """Log emergence detection events."""
        self.log_event(
            "emergence",
            universe_id,
            tick,
            {"emergence_type": emergence_type, "details": details},
            level=logging.WARNING,
        )


class PerformanceLogger:
    """Logger for performance metrics."""

    def __init__(self, name: str = "nommo.performance"):
        self.logger = logging.getLogger(name)

    def log_timing(
        self, operation: str, duration: float, extra: dict[str, Any] | None = None
    ) -> None:
        """Log operation timing."""
        data = {"operation": operation, "duration_ms": duration * 1000, **(extra or {})}
        self.logger.info(f"Performance: {operation}", extra={"extra_data": data})

    def timer(self, operation: str) -> Any:
        """Decorator for timing functions."""

        def decorator(func: Any) -> Any:
            @wraps(func)
            def wrapper(*args: Any, **kwargs: Any) -> Any:
                start = time.perf_counter()
                try:
                    result = func(*args, **kwargs)
                    duration = time.perf_counter() - start
                    self.log_timing(operation, duration, {"status": "success"})
                    return result
                except Exception as e:
                    duration = time.perf_counter() - start
                    self.log_timing(operation, duration, {"status": "error", "error": str(e)})
                    raise

            return wrapper

        return decorator


def setup_logging(
    level: str = "INFO",
    log_file: Path | None = None,
    use_json: bool = False,
    console_output: bool = True,
) -> None:
    """
    Configure logging for the entire application.

    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_file: Optional file path for file logging
        use_json: Use JSON structured logging format
        console_output: Enable console output with Rich
    """
    root_logger = logging.getLogger("nommo")
    root_logger.setLevel(getattr(logging, level.upper()))

    root_logger.handlers = []

    console_handler: logging.Handler
    if console_output and not use_json:
        console_handler = RichHandler(
            console=Console(stderr=True),
            show_time=True,
            show_path=False,
            rich_tracebacks=True,
            tracebacks_show_locals=True,
        )
        console_handler.setLevel(logging.INFO)
        root_logger.addHandler(console_handler)
    elif console_output:
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setFormatter(StructuredFormatter())
        root_logger.addHandler(console_handler)

    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.handlers.RotatingFileHandler(
            log_file, maxBytes=10 * 1024 * 1024, backupCount=5
        )

        if use_json:
            file_handler.setFormatter(StructuredFormatter())
        else:
            file_handler.setFormatter(
                logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            )

        file_handler.setLevel(logging.DEBUG)
        root_logger.addHandler(file_handler)

    logging.getLogger("numba").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("h5py").setLevel(logging.WARNING)


def get_logger(name: str) -> logging.Logger:
    """Get a logger instance for a module."""
    return logging.getLogger(f"nommo.{name}")


def get_simulation_logger() -> SimulationLogger:
    """Get the simulation event logger."""
    return SimulationLogger()


def get_performance_logger() -> PerformanceLogger:
    """Get the performance metrics logger."""
    return PerformanceLogger()
