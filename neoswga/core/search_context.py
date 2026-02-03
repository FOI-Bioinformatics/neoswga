"""
Search context dataclasses for primer optimization.

Replaces the **kwargs anti-pattern in optimize.py with strongly-typed
dataclasses that provide:
- IDE autocomplete and type checking
- Clear documentation of required vs optional parameters
- Immutable configurations that can't be accidentally modified
- Easy serialization for logging and debugging

Before:
    def bfs_one_top_set(**kwargs):
        primer_list = kwargs['primer_list']
        banned_primers = kwargs['banned_primers']
        # ... 16 more extractions with no validation

After:
    def bfs_one_top_set(ctx: BFSSearchContext) -> SearchResult:
        # All parameters typed and validated at construction
        for primer in ctx.primer_list:
            if primer not in ctx.banned_primers:
                ...
"""

from dataclasses import dataclass, field
from typing import List, Dict, Set, Tuple, Optional, Any, FrozenSet
import numpy as np
from sortedcontainers import SortedSet


@dataclass(frozen=True)
class GenomeInfo:
    """
    Immutable genome information for optimization.

    Attributes:
        prefixes: HDF5 file prefixes for this genome
        seq_lengths: Sequence lengths for each file
        circular: Whether genome is circular
    """
    prefixes: Tuple[str, ...]
    seq_lengths: Tuple[int, ...]
    circular: bool = True

    def __post_init__(self):
        if len(self.prefixes) != len(self.seq_lengths):
            raise ValueError(
                f"prefixes ({len(self.prefixes)}) and seq_lengths "
                f"({len(self.seq_lengths)}) must have same length"
            )

    @property
    def total_length(self) -> int:
        """Total genome length across all sequences."""
        return sum(self.seq_lengths)

    @property
    def num_sequences(self) -> int:
        """Number of sequences in genome."""
        return len(self.prefixes)

    @classmethod
    def from_lists(
        cls,
        prefixes: List[str],
        seq_lengths: List[int],
        circular: bool = True
    ) -> 'GenomeInfo':
        """Create from lists (converts to tuples)."""
        return cls(
            prefixes=tuple(prefixes),
            seq_lengths=tuple(seq_lengths),
            circular=circular
        )


@dataclass
class PositionData:
    """
    Mutable position data for a primer set during optimization.

    Stores binding positions for forward and reverse strands across
    all genome sequences. Updated incrementally as primers are added.
    """
    # List of (forward_positions, reverse_positions) per genome sequence
    positions: List[Tuple[SortedSet, SortedSet]]

    @classmethod
    def empty(cls, num_sequences: int) -> 'PositionData':
        """Create empty position data for given number of sequences."""
        return cls(
            positions=[(SortedSet([]), SortedSet([])) for _ in range(num_sequences)]
        )

    @classmethod
    def from_numpy(
        cls,
        position_arrays: List[Tuple[np.ndarray, np.ndarray]]
    ) -> 'PositionData':
        """Create from numpy arrays."""
        return cls(
            positions=[
                (SortedSet(fwd.tolist()), SortedSet(rev.tolist()))
                for fwd, rev in position_arrays
            ]
        )

    def copy(self) -> 'PositionData':
        """Create a deep copy."""
        return PositionData(
            positions=[
                (SortedSet(fwd), SortedSet(rev))
                for fwd, rev in self.positions
            ]
        )

    def add_positions(
        self,
        seq_idx: int,
        forward: np.ndarray,
        reverse: np.ndarray
    ) -> None:
        """Add positions for a specific sequence."""
        if seq_idx >= len(self.positions):
            raise IndexError(f"seq_idx {seq_idx} out of range")
        self.positions[seq_idx][0].update(forward.tolist())
        self.positions[seq_idx][1].update(reverse.tolist())

    def total_sites(self) -> int:
        """Total number of binding sites across all sequences."""
        return sum(
            len(fwd) + len(rev)
            for fwd, rev in self.positions
        )


@dataclass
class DimerConstraints:
    """
    Dimer checking constraints for primer compatibility.

    This class stores a pre-computed dimer compatibility matrix and provides
    efficient lookup for primer-primer interactions.

    Note: This class is mutable but should be treated as read-only after
    construction. The matrix is made read-only after initialization.

    Attributes:
        matrix: Pre-computed dimer compatibility matrix (N x N boolean)
        primer_to_index: Mapping from primer sequence to matrix index
        max_dimer_bp: Maximum allowed complementary base pairs
    """
    matrix: np.ndarray
    primer_to_index: Dict[str, int]
    max_dimer_bp: int = 4

    def __post_init__(self):
        """Make matrix read-only after construction to prevent accidental mutation."""
        if self.matrix is not None and hasattr(self.matrix, 'flags'):
            self.matrix.flags.writeable = False

    def are_compatible(self, primers: List[str], new_primer: str) -> bool:
        """
        Check if new_primer is compatible with all primers in list.

        Returns True if new_primer forms no dimers with any existing primer.
        """
        if new_primer not in self.primer_to_index:
            return True  # Unknown primer assumed compatible

        new_idx = self.primer_to_index[new_primer]

        for primer in primers:
            if primer in self.primer_to_index:
                primer_idx = self.primer_to_index[primer]
                if self.matrix[primer_idx, new_idx]:
                    return False  # Dimer detected

        return True


@dataclass
class SearchState:
    """
    Mutable state for a single search path in BFS optimization.

    Represents one of the "top sets" being tracked during search.
    """
    primers: List[str]
    score: float
    fg_positions: PositionData
    bg_positions: PositionData

    def copy(self) -> 'SearchState':
        """Create a deep copy of this state."""
        return SearchState(
            primers=list(self.primers),
            score=self.score,
            fg_positions=self.fg_positions.copy(),
            bg_positions=self.bg_positions.copy(),
        )

    def add_primer(
        self,
        primer: str,
        fg_new_positions: List[Tuple[np.ndarray, np.ndarray]],
        bg_new_positions: List[Tuple[np.ndarray, np.ndarray]],
        new_score: float
    ) -> 'SearchState':
        """Create new state with additional primer."""
        new_state = self.copy()
        new_state.primers.append(primer)
        new_state.score = new_score

        for i, (fwd, rev) in enumerate(fg_new_positions):
            new_state.fg_positions.add_positions(i, fwd, rev)

        for i, (fwd, rev) in enumerate(bg_new_positions):
            new_state.bg_positions.add_positions(i, fwd, rev)

        return new_state

    @property
    def compressed_key(self) -> str:
        """Unique string key for this primer set (for caching)."""
        return ",".join(sorted(self.primers))


@dataclass(frozen=True)
class BFSConfig:
    """
    Immutable configuration for BFS optimization.

    All parameters that control the search behavior.
    """
    max_sets: int = 10
    iterations: int = 10
    selection_method: str = 'deterministic'
    drop_indices: Tuple[int, ...] = (4,)
    normalize_metric: str = 'deterministic'
    verbose: bool = True

    def __post_init__(self):
        if self.max_sets < 1:
            raise ValueError(f"max_sets must be >= 1, got {self.max_sets}")
        if self.iterations < 1:
            raise ValueError(f"iterations must be >= 1, got {self.iterations}")
        if self.selection_method not in ('deterministic', 'softmax', 'normalized'):
            raise ValueError(f"Unknown selection_method: {self.selection_method}")


@dataclass
class BFSSearchContext:
    """
    Complete context for BFS primer set search.

    Replaces the **kwargs dictionary in bfs_one_top_set() with a
    strongly-typed, validated container.

    Usage:
        ctx = BFSSearchContext(
            primer_pool=candidates,
            fg_genome=GenomeInfo.from_lists(fg_prefixes, fg_lengths),
            bg_genome=GenomeInfo.from_lists(bg_prefixes, bg_lengths),
            dimer_constraints=DimerConstraints(matrix, primer_to_idx),
            config=BFSConfig(max_sets=10),
        )

        result = bfs_search(ctx)
    """
    # Required parameters
    primer_pool: List[str]
    fg_genome: GenomeInfo
    dimer_constraints: DimerConstraints
    config: BFSConfig

    # Optional background genome
    bg_genome: Optional[GenomeInfo] = None

    # State tracking
    banned_primers: Set[str] = field(default_factory=set)
    score_cache: Dict[str, float] = field(default_factory=dict)

    # Current search state (mutable during search)
    top_states: List[SearchState] = field(default_factory=list)

    def __post_init__(self):
        # Validate primer pool
        if not self.primer_pool:
            raise ValueError("primer_pool cannot be empty")

        # Deduplicate primer pool
        seen = set()
        unique = []
        for p in self.primer_pool:
            if p not in seen:
                seen.add(p)
                unique.append(p)
        self.primer_pool = unique

    @property
    def has_background(self) -> bool:
        """Whether background genome is configured."""
        return self.bg_genome is not None

    def get_available_primers(self, current_set: List[str]) -> List[str]:
        """Get primers that can be added to current set."""
        current_set_frozen = frozenset(current_set)
        return [
            p for p in self.primer_pool
            if p not in self.banned_primers
            and p not in current_set_frozen
            and self.dimer_constraints.are_compatible(current_set, p)
        ]

    def get_cached_score(self, primers: List[str]) -> Optional[float]:
        """Get cached score for primer set if available."""
        key = ",".join(sorted(primers))
        return self.score_cache.get(key)

    def cache_score(self, primers: List[str], score: float) -> None:
        """Cache score for primer set."""
        key = ",".join(sorted(primers))
        self.score_cache[key] = score

    def initialize_states(
        self,
        initial_sets: Optional[List[List[str]]] = None
    ) -> None:
        """
        Initialize search states.

        Args:
            initial_sets: Optional initial primer sets to start from.
                         If None, starts with empty sets.
        """
        if initial_sets:
            # Build states from initial sets
            self.top_states = []
            for primers in initial_sets[:self.config.max_sets]:
                state = SearchState(
                    primers=list(primers),
                    score=float('-inf'),
                    fg_positions=PositionData.empty(self.fg_genome.num_sequences),
                    bg_positions=PositionData.empty(
                        self.bg_genome.num_sequences if self.bg_genome else 0
                    ),
                )
                self.top_states.append(state)
        else:
            # Start with empty sets
            self.top_states = [
                SearchState(
                    primers=[],
                    score=float('-inf'),
                    fg_positions=PositionData.empty(self.fg_genome.num_sequences),
                    bg_positions=PositionData.empty(
                        self.bg_genome.num_sequences if self.bg_genome else 0
                    ),
                )
                for _ in range(self.config.max_sets)
            ]


@dataclass
class SearchResult:
    """
    Result from a single search iteration or complete search.
    """
    states: List[SearchState]
    best_score: float
    iterations_completed: int
    converged: bool = False
    message: str = ""

    @property
    def best_state(self) -> Optional[SearchState]:
        """Get the highest-scoring state."""
        if not self.states:
            return None
        return max(self.states, key=lambda s: s.score)

    @property
    def best_primers(self) -> List[str]:
        """Get primers from the best state."""
        best = self.best_state
        return best.primers if best else []

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'best_primers': self.best_primers,
            'best_score': self.best_score,
            'iterations': self.iterations_completed,
            'converged': self.converged,
            'message': self.message,
            'num_states': len(self.states),
        }


@dataclass(frozen=True)
class EvaluationContext:
    """
    Immutable context for primer set evaluation.

    Used by the evaluate() function to score primer sets.
    """
    fg_genome: GenomeInfo
    bg_genome: Optional[GenomeInfo] = None

    # Model path for ridge regression
    model_path: Optional[str] = None

    def has_background(self) -> bool:
        """Whether background genome is available."""
        return self.bg_genome is not None
