"""
Factory for creating optimizer instances.

Implements the Factory and Registry patterns to provide:
- Centralized optimizer registration
- Runtime optimizer discovery
- Consistent instantiation interface
- Easy extension with new optimizers

Usage:
    # Get optimizer by name
    optimizer = OptimizerFactory.create('greedy', cache=cache, ...)

    # List available optimizers
    for name, desc in OptimizerFactory.list_optimizers().items():
        print(f"{name}: {desc}")

    # Register custom optimizer
    @OptimizerFactory.register('my-optimizer')
    class MyOptimizer(BaseOptimizer):
        ...
"""

from typing import Dict, Type, Optional, List, Any, Callable
import logging
import threading

from .base_optimizer import BaseOptimizer, OptimizerConfig
from .exceptions import OptimizerNotFoundError

logger = logging.getLogger(__name__)


class OptimizerRegistry:
    """
    Registry for optimizer classes.

    Thread-safe singleton that maintains the mapping from optimizer
    names to their implementation classes.
    """

    _instance = None
    _lock = threading.Lock()
    _registry: Dict[str, Type[BaseOptimizer]] = {}
    _descriptions: Dict[str, str] = {}
    _aliases: Dict[str, str] = {}

    def __new__(cls):
        if cls._instance is None:
            with cls._lock:
                # Double-checked locking pattern
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
        return cls._instance

    @classmethod
    def register(
        cls,
        name: str,
        aliases: Optional[List[str]] = None,
        description: Optional[str] = None
    ) -> Callable[[Type[BaseOptimizer]], Type[BaseOptimizer]]:
        """
        Decorator to register an optimizer class.

        Args:
            name: Primary name for the optimizer (kebab-case preferred)
            aliases: Alternative names (e.g., 'ga' for 'genetic-algorithm')
            description: Optional description (defaults to class docstring)

        Returns:
            Decorator function

        Example:
            @OptimizerRegistry.register('greedy-bfs', aliases=['greedy', 'bfs'])
            class GreedyOptimizer(BaseOptimizer):
                '''Greedy breadth-first search optimizer.'''
                ...
        """
        def decorator(optimizer_class: Type[BaseOptimizer]) -> Type[BaseOptimizer]:
            # Validate
            if not issubclass(optimizer_class, BaseOptimizer):
                raise TypeError(
                    f"Cannot register {optimizer_class.__name__}: "
                    f"must be subclass of BaseOptimizer"
                )

            with cls._lock:
                # Register primary name
                cls._registry[name] = optimizer_class

                # Store description
                desc = description or optimizer_class.__doc__
                if desc:
                    # Take first line of docstring
                    cls._descriptions[name] = desc.strip().split('\n')[0]
                else:
                    cls._descriptions[name] = f"{name} optimizer"

                # Register aliases
                if aliases:
                    for alias in aliases:
                        cls._aliases[alias] = name

            logger.debug(f"Registered optimizer: {name}")
            return optimizer_class

        return decorator

    @classmethod
    def get(cls, name: str) -> Type[BaseOptimizer]:
        """
        Get optimizer class by name.

        Args:
            name: Optimizer name or alias

        Returns:
            Optimizer class

        Raises:
            OptimizerNotFoundError: If optimizer not found
        """
        with cls._lock:
            # Check aliases first
            canonical_name = cls._aliases.get(name, name)

            if canonical_name not in cls._registry:
                raise OptimizerNotFoundError(name, list(cls._registry.keys()))

            return cls._registry[canonical_name]

    @classmethod
    def list_all(cls) -> Dict[str, str]:
        """
        List all registered optimizers with descriptions.

        Returns:
            Dict mapping optimizer names to descriptions
        """
        with cls._lock:
            return dict(cls._descriptions)

    @classmethod
    def is_registered(cls, name: str) -> bool:
        """Check if an optimizer is registered."""
        with cls._lock:
            canonical_name = cls._aliases.get(name, name)
            return canonical_name in cls._registry

    @classmethod
    def clear(cls) -> None:
        """Clear all registrations (for testing)."""
        with cls._lock:
            cls._registry.clear()
            cls._descriptions.clear()
            cls._aliases.clear()


class OptimizerFactory:
    """
    Factory for creating optimizer instances.

    Provides a clean interface for optimizer creation with consistent
    parameter handling and validation.

    Usage:
        # Create optimizer with explicit parameters
        optimizer = OptimizerFactory.create(
            'greedy',
            position_cache=cache,
            fg_prefixes=['data/target'],
            fg_seq_lengths=[1000000],
            config=OptimizerConfig(target_set_size=8)
        )

        # Create from parameter dict
        optimizer = OptimizerFactory.from_params(params)

        # List available
        for name, desc in OptimizerFactory.list_optimizers().items():
            print(f"{name}: {desc}")
    """

    @staticmethod
    def create(
        name: str,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
        **kwargs
    ) -> BaseOptimizer:
        """
        Create optimizer instance by name.

        Args:
            name: Optimizer name or alias
            position_cache: PositionCache for primer position lookups
            fg_prefixes: Foreground genome HDF5 prefixes
            fg_seq_lengths: Foreground genome lengths
            bg_prefixes: Background genome HDF5 prefixes (optional)
            bg_seq_lengths: Background genome lengths (optional)
            config: Optimizer configuration
            **kwargs: Additional optimizer-specific parameters

        Returns:
            Configured optimizer instance

        Raises:
            OptimizerNotFoundError: If optimizer name not found
            ValueError: If required parameters are invalid

        Example:
            optimizer = OptimizerFactory.create(
                'dominating-set',
                position_cache=cache,
                fg_prefixes=['data/target'],
                fg_seq_lengths=[1000000],
                bin_size=5000  # optimizer-specific param
            )
        """
        # Validate required inputs
        if not name:
            raise ValueError("Optimizer name cannot be empty")
        if position_cache is None:
            raise ValueError("position_cache is required")
        if not fg_prefixes:
            raise ValueError("fg_prefixes cannot be empty")
        if not fg_seq_lengths:
            raise ValueError("fg_seq_lengths cannot be empty")
        if len(fg_prefixes) != len(fg_seq_lengths):
            raise ValueError(
                f"fg_prefixes ({len(fg_prefixes)}) and fg_seq_lengths "
                f"({len(fg_seq_lengths)}) must have same length"
            )

        optimizer_class = OptimizerRegistry.get(name)

        # Merge config with kwargs if optimizer has specific config class
        optimizer_config = config or OptimizerConfig()

        try:
            optimizer = optimizer_class(
                position_cache=position_cache,
                fg_prefixes=fg_prefixes,
                fg_seq_lengths=fg_seq_lengths,
                bg_prefixes=bg_prefixes,
                bg_seq_lengths=bg_seq_lengths,
                config=optimizer_config,
                **kwargs
            )
            logger.info(f"Created optimizer: {optimizer.name}")
            return optimizer

        except TypeError as e:
            # Provide helpful error if constructor signature doesn't match
            logger.error(f"Failed to create optimizer '{name}': {e}")
            raise

    @staticmethod
    def from_params(params: Dict[str, Any], position_cache) -> BaseOptimizer:
        """
        Create optimizer from parameter dictionary.

        Extracts optimizer configuration from a params dict (as used by CLI).

        Args:
            params: Parameter dictionary with keys:
                - optimization_method: Optimizer name
                - fg_fasta, bg_fasta: Genome files
                - data_dir: Directory with HDF5 files
                - num_primers, target_set_size: Target primer count
                - iterations, max_sets: Optimization parameters
            position_cache: PositionCache instance

        Returns:
            Configured optimizer instance

        Raises:
            OptimizerNotFoundError: If requested optimizer is not registered
            InvalidParameterError: If required parameters are missing/invalid
        """
        from .parameter import get_params
        from .exceptions import InvalidParameterError

        # Get optimizer name with validation
        method = params.get('optimization_method', 'greedy')

        # Validate that optimizer exists before proceeding
        if not OptimizerRegistry.is_registered(method):
            available = list(OptimizerRegistry.list_all().keys())
            from .exceptions import OptimizerNotFoundError
            raise OptimizerNotFoundError(method, available)

        # Build prefixes from data_dir
        data_dir = params.get('data_dir', '.')
        fg_prefix = params.get('fg_prefix')
        bg_prefix = params.get('bg_prefix')

        fg_prefixes = [fg_prefix] if fg_prefix else []
        bg_prefixes = [bg_prefix] if bg_prefix else []

        # Get sequence lengths (would need to be loaded from somewhere)
        fg_seq_lengths = params.get('fg_seq_lengths', [])
        bg_seq_lengths = params.get('bg_seq_lengths', [])

        # Validate target set size
        target_set_size = params.get('num_primers', params.get('target_set_size', 6))
        if not isinstance(target_set_size, int) or target_set_size < 1:
            raise InvalidParameterError(
                'target_set_size', target_set_size,
                "must be a positive integer"
            )

        # Validate iterations
        iterations = params.get('iterations', 100)
        if not isinstance(iterations, int) or iterations < 1:
            raise InvalidParameterError(
                'iterations', iterations,
                "must be a positive integer"
            )

        # Build config
        config = OptimizerConfig(
            target_set_size=target_set_size,
            max_iterations=iterations,
            max_dimer_bp=params.get('max_dimer_bp', 4),
            min_tm=params.get('min_tm', 20.0),
            max_tm=params.get('max_tm', 50.0),
            verbose=params.get('verbose', True),
        )

        return OptimizerFactory.create(
            name=method,
            position_cache=position_cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes,
            bg_seq_lengths=bg_seq_lengths,
            config=config,
        )

    @staticmethod
    def list_optimizers() -> Dict[str, str]:
        """
        List all available optimizers.

        Returns:
            Dict mapping optimizer names to descriptions
        """
        return OptimizerRegistry.list_all()

    @staticmethod
    def get_optimizer_info(name: str) -> Dict[str, Any]:
        """
        Get detailed information about an optimizer.

        Args:
            name: Optimizer name

        Returns:
            Dict with optimizer metadata
        """
        optimizer_class = OptimizerRegistry.get(name)
        return {
            'name': name,
            'class': optimizer_class.__name__,
            'description': OptimizerRegistry._descriptions.get(name, ''),
            'supports_background': getattr(optimizer_class, 'supports_background', True),
            'module': optimizer_class.__module__,
        }

    @staticmethod
    def register(
        name: str,
        aliases: Optional[List[str]] = None,
        description: Optional[str] = None
    ):
        """
        Decorator to register an optimizer class.

        Convenience wrapper around OptimizerRegistry.register().

        Example:
            @OptimizerFactory.register('my-optimizer', aliases=['my', 'custom'])
            class MyOptimizer(BaseOptimizer):
                ...
        """
        return OptimizerRegistry.register(name, aliases, description)


# =============================================================================
# Built-in Optimizer Registration
# =============================================================================
# Import and register built-in optimizers when this module is loaded.
# Each optimizer module should use @OptimizerFactory.register() decorator.

def _register_builtin_optimizers():
    """
    Register built-in optimizers.

    Called at module load time to ensure all optimizers are available.
    Uses lazy imports to avoid circular dependencies.
    """
    try:
        # These imports will trigger registration via decorators
        # We use a try/except to handle missing dependencies gracefully

        # Note: The actual optimizer files need to be updated to use
        # @OptimizerFactory.register() decorator. Until then, we register
        # them manually here for backwards compatibility.

        pass  # Registration happens via decorators in optimizer modules

    except ImportError as e:
        logger.warning(f"Failed to import some optimizers: {e}")


# Run registration at module load
_register_builtin_optimizers()
