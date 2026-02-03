"""
Progress indicators for neoswga pipeline.

Provides simple, terminal-friendly progress display for long-running operations.
Works without external dependencies (no tqdm required).
"""

import sys
import time
from typing import Optional, Iterable, TypeVar, Iterator
from contextlib import contextmanager

T = TypeVar('T')


class ProgressBar:
    """
    Simple progress bar for terminal output.

    Usage:
        # With known total
        with ProgressBar(total=100, desc="Processing") as pbar:
            for item in items:
                process(item)
                pbar.update(1)

        # With iterable
        for item in ProgressBar.wrap(items, desc="Processing"):
            process(item)
    """

    def __init__(
        self,
        total: Optional[int] = None,
        desc: str = "",
        unit: str = "it",
        width: int = 40,
        disable: bool = False
    ):
        """
        Initialize progress bar.

        Args:
            total: Total number of items (None for unknown)
            desc: Description text
            unit: Unit name for items
            width: Bar width in characters
            disable: Disable progress output
        """
        self.total = total
        self.desc = desc
        self.unit = unit
        self.width = width
        self.disable = disable

        self.n = 0
        self.start_time = None
        self.last_print_time = 0

    def __enter__(self):
        self.start_time = time.time()
        self._print()
        return self

    def __exit__(self, *args):
        self._print(force=True, newline=True)

    def update(self, n: int = 1) -> None:
        """Update progress by n items."""
        self.n += n

        # Rate-limit updates to avoid flickering
        now = time.time()
        if now - self.last_print_time >= 0.1:
            self._print()
            self.last_print_time = now

    def set_description(self, desc: str) -> None:
        """Update description text."""
        self.desc = desc

    def _print(self, force: bool = False, newline: bool = False) -> None:
        """Print progress bar."""
        if self.disable:
            return

        elapsed = time.time() - self.start_time if self.start_time else 0

        if self.total:
            # Known total - show percentage bar
            pct = self.n / self.total
            filled = int(self.width * pct)
            bar = "=" * filled + ">" + " " * (self.width - filled - 1)

            rate = self.n / elapsed if elapsed > 0 else 0
            eta = (self.total - self.n) / rate if rate > 0 else 0

            line = f"\r{self.desc}: [{bar}] {self.n}/{self.total} ({pct:.0%}) | {rate:.1f} {self.unit}/s | ETA: {_format_time(eta)}"
        else:
            # Unknown total - show spinner and count
            spinner = "|/-\\"[int(elapsed * 4) % 4]
            rate = self.n / elapsed if elapsed > 0 else 0
            line = f"\r{self.desc}: {spinner} {self.n} {self.unit} | {rate:.1f} {self.unit}/s | {_format_time(elapsed)}"

        sys.stdout.write(line)
        sys.stdout.flush()

        if newline:
            sys.stdout.write("\n")

    @classmethod
    def wrap(
        cls,
        iterable: Iterable[T],
        total: Optional[int] = None,
        desc: str = "",
        unit: str = "it",
        disable: bool = False
    ) -> Iterator[T]:
        """
        Wrap an iterable with a progress bar.

        Args:
            iterable: Iterable to wrap
            total: Total count (auto-detected if possible)
            desc: Description text
            unit: Unit name
            disable: Disable progress

        Yields:
            Items from the iterable
        """
        if total is None:
            try:
                total = len(iterable)
            except TypeError:
                pass

        with cls(total=total, desc=desc, unit=unit, disable=disable) as pbar:
            for item in iterable:
                yield item
                pbar.update(1)


class StepProgress:
    """
    Progress indicator for multi-step operations.

    Usage:
        progress = StepProgress(steps=[
            "Loading genome",
            "Counting k-mers",
            "Writing output"
        ])

        progress.start_step(0)
        # ... do work ...
        progress.complete_step(0)

        progress.start_step(1)
        # ... do work ...
        progress.complete_step(1)
    """

    def __init__(self, steps: list, disable: bool = False):
        """
        Initialize step progress.

        Args:
            steps: List of step descriptions
            disable: Disable output
        """
        self.steps = steps
        self.disable = disable
        self.current_step = -1
        self.step_status = ['pending'] * len(steps)
        self.start_time = time.time()

    def start_step(self, step: int) -> None:
        """Mark a step as started."""
        self.current_step = step
        self.step_status[step] = 'running'
        if not self.disable:
            self._print_status()

    def complete_step(self, step: int) -> None:
        """Mark a step as completed."""
        self.step_status[step] = 'done'
        if not self.disable:
            self._print_status()

    def fail_step(self, step: int, error: str = "") -> None:
        """Mark a step as failed."""
        self.step_status[step] = 'failed'
        if not self.disable:
            self._print_status()
            if error:
                print(f"  Error: {error}")

    def _print_status(self) -> None:
        """Print current status."""
        elapsed = time.time() - self.start_time

        print(f"\n--- Pipeline Progress ({_format_time(elapsed)}) ---")
        for i, (step, status) in enumerate(zip(self.steps, self.step_status)):
            if status == 'done':
                icon = "[OK]"
            elif status == 'running':
                icon = "[..]"
            elif status == 'failed':
                icon = "[!!]"
            else:
                icon = "[  ]"

            print(f"  {icon} Step {i+1}: {step}")


@contextmanager
def progress_context(desc: str, disable: bool = False):
    """
    Context manager for simple start/complete progress.

    Usage:
        with progress_context("Loading genome"):
            load_genome()
    """
    start = time.time()  # Always capture start time for safety
    if not disable:
        print(f"{desc}...", end=" ", flush=True)

    try:
        yield
        if not disable:
            elapsed = time.time() - start
            print(f"done ({elapsed:.1f}s)")
    except Exception:
        if not disable:
            print("FAILED")
        raise


def _format_time(seconds: float) -> str:
    """Format seconds as human-readable time."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        m, s = divmod(int(seconds), 60)
        return f"{m}m {s}s"
    else:
        h, rem = divmod(int(seconds), 3600)
        m, s = divmod(rem, 60)
        return f"{h}h {m}m"


# Convenience functions
def progress_bar(iterable, **kwargs):
    """Shortcut for ProgressBar.wrap()."""
    return ProgressBar.wrap(iterable, **kwargs)


def step_progress(steps, **kwargs):
    """Shortcut for StepProgress()."""
    return StepProgress(steps, **kwargs)
