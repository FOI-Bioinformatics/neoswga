"""
Network-based primer set optimization.

KEY INSIGHT: SWGA amplification efficiency depends on NETWORK CONNECTIVITY,
not just binding frequency.

Connected primers → Exponential amplification
Isolated primers → Linear amplification

This is the fundamental difference from the ratio-based approach.
"""

import math
import numpy as np
import networkx as nx
from typing import List, Dict, Tuple, Set, Optional, TYPE_CHECKING
from dataclasses import dataclass, field
import logging

if TYPE_CHECKING:
    from neoswga.core.reaction_conditions import ReactionConditions
from collections import defaultdict
import bisect

logger = logging.getLogger(__name__)


def calculate_primer_tm(primer: str) -> float:
    """
    Calculate melting temperature using nearest-neighbor thermodynamics.
    Falls back to simple formula if melting module unavailable.
    """
    from neoswga.core.melting_temp import temp as _melting_temp
    return _melting_temp(primer)


def calculate_dimer_score(primer1: str, primer2: str, max_bp: int = 4) -> float:
    """
    Calculate dimer formation potential between two primers.
    Returns score 0-1 where 0 = no dimer risk, 1 = high dimer risk.
    """
    try:
        from neoswga.core import dimer
        if dimer.is_dimer_fast(primer1, primer2, max_bp):
            return 1.0
        return 0.0
    except ImportError:
        # Simplified dimer check: look for complementary runs
        _comp_table = str.maketrans('ATGC', 'TACG')
        p1 = primer1.upper()
        p2_rev = primer2.upper()[::-1]
        p2_comp = p2_rev.translate(_comp_table)

        # Check for runs of complementary bases
        max_run = 0
        for i in range(len(p1)):
            for j in range(len(p2_comp)):
                run = 0
                while i + run < len(p1) and j + run < len(p2_comp):
                    if p1[i + run] == p2_comp[j + run]:
                        run += 1
                    else:
                        break
                max_run = max(max_run, run)

        return min(1.0, max_run / max_bp) if max_run >= max_bp else 0.0


@dataclass(slots=True)
class BindingSite:
    """Primer binding site"""
    position: int
    strand: str  # '+' or '-'
    primer: str
    affinity: float = 1.0  # Binding affinity (1.0 = perfect match)
    _hash: int = field(default=0, init=False, repr=False, compare=False)

    def __hash__(self):
        h = self._hash
        if h == 0:
            h = hash((self.position, self.strand, self.primer))
            self._hash = h
        return h


@dataclass
class AmplificationEdge:
    """Edge in amplification network (extension path)"""
    source: BindingSite
    target: BindingSite
    distance: int
    feasible: bool  # Can Phi29 extend this far?

    def __hash__(self):
        return hash((self.source, self.target))


class AmplificationNetwork:
    """
    Graph representation of primer amplification potential.

    Nodes: Binding sites (primer locations)
    Edges: Extension paths (Phi29 can reach from site A to site B)

    Network structure determines amplification:
    - Large connected component → Exponential growth
    - Small disconnected components → Linear growth

    Uses spatial indexing for efficient edge building (O(n log n) instead of O(n^2)).
    """

    def __init__(self, max_extension: int = 70000):
        """
        Initialize network.

        Args:
            max_extension: Maximum Phi29 extension length (default 70 kb)
        """
        self.graph = nx.Graph()
        self.max_extension = max_extension
        self.binding_sites: List[BindingSite] = []

        # Spatial index: sorted positions and corresponding sites for binary search
        self._sorted_positions: List[int] = []
        self._sorted_sites: List[BindingSite] = []
        self._index_dirty: bool = True

    def _rebuild_spatial_index(self):
        """Rebuild sorted spatial index for efficient range queries."""
        if not self._index_dirty:
            return

        # Sort sites by position
        sorted_pairs = sorted(
            [(s.position, s) for s in self.binding_sites],
            key=lambda x: x[0]
        )
        self._sorted_positions = [p for p, _ in sorted_pairs]
        self._sorted_sites = [s for _, s in sorted_pairs]
        self._index_dirty = False

    def _find_sites_in_range(self, pos: int, max_dist: int) -> List[BindingSite]:
        """
        Find all sites within max_dist of position using binary search.

        Returns sites in range [pos - max_dist, pos + max_dist].
        O(log n + k) where k is number of sites in range.
        """
        self._rebuild_spatial_index()

        if not self._sorted_positions:
            return []

        # Find range bounds using bisect
        left_bound = pos - max_dist
        right_bound = pos + max_dist

        left_idx = bisect.bisect_left(self._sorted_positions, left_bound)
        right_idx = bisect.bisect_right(self._sorted_positions, right_bound)

        return self._sorted_sites[left_idx:right_idx]

    def add_primer_sites(self, primer: str, positions: np.ndarray,
                        strand: str, affinity: float = 1.0) -> List[BindingSite]:
        """
        Add all binding sites for a primer.

        Args:
            primer: Primer sequence
            positions: Array of binding positions
            strand: '+' or '-'
            affinity: Binding strength (1.0 = perfect)

        Returns:
            List of newly added BindingSite objects
        """
        new_sites = []
        for pos in positions:
            site = BindingSite(position=int(pos), strand=strand,
                             primer=primer, affinity=affinity)
            self.binding_sites.append(site)
            self.graph.add_node(site)
            new_sites.append(site)

        # Mark spatial index as dirty
        self._index_dirty = True
        return new_sites

    def add_edges_for_sites(self, new_sites: List[BindingSite]) -> int:
        """
        Add edges only for newly added sites (incremental edge building).

        Much faster than rebuild_edges() when adding a few sites to large network.
        O(k * log n + k * m) where k = new sites, m = avg sites in range.

        Args:
            new_sites: List of newly added binding sites

        Returns:
            Number of edges added
        """
        self._rebuild_spatial_index()
        edges_added = 0

        for new_site in new_sites:
            # Find all existing sites within extension distance
            nearby_sites = self._find_sites_in_range(
                new_site.position, self.max_extension
            )

            for other_site in nearby_sites:
                # Skip self
                if other_site is new_site:
                    continue

                # Must be opposite strands
                if new_site.strand == other_site.strand:
                    continue

                # Check if edge already exists
                if self.graph.has_edge(new_site, other_site):
                    continue

                distance = abs(new_site.position - other_site.position)
                edge = AmplificationEdge(
                    source=new_site,
                    target=other_site,
                    distance=distance,
                    feasible=True
                )
                self.graph.add_edge(new_site, other_site, edge=edge, distance=distance)
                edges_added += 1

        return edges_added

    def build_edges(self):
        """
        Build edges based on extension feasibility.

        Two sites connect if:
        1. On opposite strands (required for extension)
        2. Within max_extension distance

        Uses spatial indexing for O(n log n) performance.
        """
        n_sites = len(self.binding_sites)
        if n_sites > 1000:
            logger.info(f"Building amplification network: {n_sites} sites")

        self._rebuild_spatial_index()
        edges_added = 0

        # Use sliding window on sorted sites
        for i, site1 in enumerate(self._sorted_sites):
            # Binary search for end of valid range
            max_pos = site1.position + self.max_extension
            end_idx = bisect.bisect_right(self._sorted_positions, max_pos)

            # Check sites in range
            for j in range(i + 1, end_idx):
                site2 = self._sorted_sites[j]

                # Must be opposite strands for amplification
                if site1.strand != site2.strand:
                    distance = site2.position - site1.position
                    edge = AmplificationEdge(
                        source=site1,
                        target=site2,
                        distance=distance,
                        feasible=True
                    )
                    self.graph.add_edge(site1, site2, edge=edge, distance=distance)
                    edges_added += 1

        if n_sites > 1000:
            logger.info(f"Added {edges_added} amplification edges")

    def largest_component_size(self) -> int:
        """Size of largest connected component"""
        if len(self.graph.nodes()) == 0:
            return 0
        components = nx.connected_components(self.graph)
        return max(len(c) for c in components)

    def average_component_size(self) -> float:
        """Average size of connected components"""
        if len(self.graph.nodes()) == 0:
            return 0.0
        components = list(nx.connected_components(self.graph))
        return np.mean([len(c) for c in components])

    def num_components(self) -> int:
        """Number of connected components"""
        if len(self.graph.nodes()) == 0:
            return 0
        return nx.number_connected_components(self.graph)

    def connectivity_score(self) -> float:
        """
        Overall network connectivity (algebraic connectivity).

        Higher value = better connected network.
        Zero = disconnected graph.
        """
        if len(self.graph.nodes()) < 2:
            return 0.0

        try:
            # Fiedler value (second-smallest eigenvalue of Laplacian).
            # Use tracemin_lu (direct solver) because the default tracemin_pcg
            # fails to converge on bipartite-like amplification networks.
            return nx.algebraic_connectivity(self.graph, method='tracemin_lu')
        except (nx.NetworkXError, ValueError, np.linalg.LinAlgError):
            # Graph not connected or numerical issues
            return 0.0

    def coverage_uniformity(self, genome_length: int) -> float:
        """
        Measure how uniformly network covers genome.

        Returns CV (coefficient of variation) of component positions.
        Lower is better.
        """
        components = list(nx.connected_components(self.graph))

        if len(components) == 0:
            return float('inf')

        # Average position of each component
        component_positions = []
        for component in components:
            positions = [site.position for site in component]
            component_positions.append(np.mean(positions))

        if len(component_positions) < 2:
            return 0.0

        # CV of component positions (coefficient of variation)
        mean_pos = np.mean(component_positions)
        if mean_pos == 0:
            return 0.0  # Avoid division by zero
        return np.std(component_positions) / mean_pos

    def predict_amplification_fold(self) -> float:
        """
        Predict amplification fold based on network structure.

        Exponential growth: 2^(component_size / 10)
        Cap at 2^20 (realistic maximum)
        """
        largest = self.largest_component_size()

        if largest < 10:
            # Too small for exponential
            return largest * 5  # Linear growth

        # Exponential scaling with saturation
        exponent = min(largest / 10.0, 20.0)
        return 2 ** exponent

    def predict_amplification_fold_weighted(self, primer_tms: Dict[str, float],
                                           reaction_temp: float) -> float:
        """
        Predict amplification accounting for Tm efficiency.

        Primers with Tm far from reaction temp contribute less to effective
        network connectivity due to reduced binding efficiency.

        Args:
            primer_tms: Dictionary mapping primer sequence to its melting temperature
            reaction_temp: Reaction temperature in Celsius

        Returns:
            Predicted amplification fold with Tm weighting
        """
        if not self.binding_sites:
            return 1.0

        # Calculate efficiency for each primer based on Tm vs reaction temp
        efficiencies = {}
        for primer, tm in primer_tms.items():
            # Bell curve: efficiency drops as |Tm - reaction_temp| increases
            # At optimal Tm (reaction_temp + 5-10C), efficiency is ~1.0
            # Efficiency drops to ~50% at 3C difference
            delta = abs(tm - (reaction_temp + 5))  # Optimal is ~5C above reaction temp
            efficiency = math.exp(-0.1 * delta**2)  # Gaussian decay
            efficiencies[primer] = max(0.1, efficiency)  # Floor at 10%

        # Weight connected component by average primer efficiency
        components = list(nx.connected_components(self.graph))
        if not components:
            return 1.0

        largest_component = max(components, key=len)
        weighted_size = 0.0

        for site in largest_component:
            eff = efficiencies.get(site.primer, 0.5)
            weighted_size += eff * site.affinity

        # Apply exponential model with weighted size
        if weighted_size < 10:
            return weighted_size * 5  # Linear for small networks

        exponent = min(weighted_size / 10.0, 20.0)
        return 2 ** exponent

    def predict_amplification_with_uncertainty(self) -> Tuple[float, float, float]:
        """
        Return (predicted, lower_bound, upper_bound) for amplification.

        Uncertainty increases with:
        - Fewer binding sites (statistical uncertainty)
        - Higher Gini index (uneven coverage)
        - Smaller network (less robust)

        Returns:
            Tuple of (point_estimate, lower_95_ci, upper_95_ci)
        """
        base_pred = self.predict_amplification_fold()

        if base_pred <= 1.0:
            return (1.0, 1.0, 1.0)

        # Calculate uncertainty factors
        num_sites = len(self.binding_sites)

        # Statistical uncertainty: fewer sites = more uncertainty
        # sqrt(N) scaling for standard error
        site_factor = 1.0 / math.sqrt(max(num_sites, 1))

        # Coverage uniformity factor
        # Higher Gini = more uneven = more uncertainty
        gini = self._calculate_gini_coefficient()
        gini_factor = gini * 0.5  # Scale Gini contribution

        # Network structure factor
        # Smaller components = less robust amplification
        largest = self.largest_component_size()
        structure_factor = max(0, 1.0 - largest / 100.0)  # Less uncertain with larger components

        # Combined uncertainty (multiplicative in log space)
        # Base uncertainty of 20% plus contributions from factors
        log_uncertainty = 0.2 + site_factor + gini_factor + structure_factor * 0.3

        # Convert to confidence interval (assuming log-normal distribution)
        log_pred = math.log(base_pred)
        lower = math.exp(log_pred - 1.96 * log_uncertainty * log_pred)
        upper = math.exp(log_pred + 1.96 * log_uncertainty * log_pred)

        # Ensure sensible bounds
        lower = max(1.0, lower)
        upper = min(2**20, upper)  # Cap at realistic maximum

        return (base_pred, lower, upper)

    def _calculate_gini_coefficient(self) -> float:
        """
        Calculate Gini coefficient of binding site positions.

        Returns 0-1 where 0 = perfectly uniform, 1 = maximally uneven.

        Uses vectorized numpy operations for O(n log n) instead of O(n^2).
        """
        if len(self.binding_sites) < 2:
            return 0.0

        positions = np.array([s.position for s in self.binding_sites])
        n = len(positions)

        # Sort positions
        positions = np.sort(positions)

        # Mean absolute difference formula for Gini (vectorized)
        # Gini = (2 * sum(i * x_i) - (n + 1) * sum(x_i)) / (n * sum(x_i))
        # For sorted data: Gini = (2 * sum((i+1) * x_i)) / (n * sum(x_i)) - (n+1)/n
        total = positions.sum()
        if total == 0:
            return 0.0

        # Weighted sum where weights are positions in sorted order (1, 2, 3, ...)
        weighted_sum = np.sum((np.arange(1, n + 1)) * positions)

        gini = (2.0 * weighted_sum) / (n * total) - (n + 1.0) / n

        return min(1.0, max(0.0, gini))

    def get_statistics(self) -> Dict:
        """Get comprehensive network statistics including confidence intervals.

        Computes connected components once and reuses across all metrics.
        """
        if len(self.graph.nodes()) == 0:
            return {
                'num_sites': 0,
                'num_edges': 0,
                'num_components': 0,
                'largest_component': 0,
                'avg_component_size': 0,
                'connectivity': 0,
                'predicted_amplification': 0,
                'amplification_lower_ci': 0,
                'amplification_upper_ci': 0,
                'gini_coefficient': 0,
            }

        # Compute components once, reuse for all derived metrics
        components = list(nx.connected_components(self.graph))
        comp_sizes = [len(c) for c in components]
        largest = max(comp_sizes)
        avg_size = np.mean(comp_sizes)
        gini = self._calculate_gini_coefficient()

        # Inline amplification prediction to avoid redundant component calls
        if largest < 10:
            base_pred = largest * 5
        else:
            exponent = min(largest / 10.0, 20.0)
            base_pred = 2 ** exponent

        if base_pred <= 1.0:
            pred, lower, upper = 1.0, 1.0, 1.0
        else:
            num_sites = len(self.binding_sites)
            site_factor = 1.0 / math.sqrt(max(num_sites, 1))
            gini_factor = gini * 0.5
            structure_factor = max(0, 1.0 - largest / 100.0)
            log_uncertainty = 0.2 + site_factor + gini_factor + structure_factor * 0.3
            log_pred = math.log(base_pred)
            lower = max(1.0, math.exp(log_pred - 1.96 * log_uncertainty * log_pred))
            upper = min(2**20, math.exp(log_pred + 1.96 * log_uncertainty * log_pred))
            pred = base_pred

        return {
            'num_sites': len(self.graph.nodes()),
            'num_edges': len(self.graph.edges()),
            'num_components': len(components),
            'largest_component': largest,
            'avg_component_size': avg_size,
            'connectivity': self.connectivity_score(),
            'predicted_amplification': pred,
            'amplification_lower_ci': lower,
            'amplification_upper_ci': upper,
            'gini_coefficient': gini,
        }


class _SimulatedNetwork:
    """Lightweight proxy storing pre-computed metrics from _simulate_add_primer.

    Provides the same interface as AmplificationNetwork for the subset of
    methods used by _evaluate_primer_addition, but without holding a full
    graph copy.
    """

    __slots__ = ('_largest', '_avg', '_uniformity', '_has_sites')

    def __init__(self, largest: int, avg: float, uniformity: float, has_sites: bool):
        self._largest = largest
        self._avg = avg
        self._uniformity = uniformity
        self._has_sites = has_sites

    def largest_component_size(self) -> int:
        return self._largest

    def average_component_size(self) -> float:
        return self._avg

    def coverage_uniformity(self, genome_length: int) -> float:
        return self._uniformity

    @property
    def binding_sites(self):
        """Return truthy/falsy value matching the has-sites check in callers."""
        return self._has_sites


class NetworkOptimizer:
    """
    Optimize primer sets for network connectivity.

    This is the core algorithm replacing greedy BFS.
    """

    def __init__(self, position_cache, fg_prefixes: List[str],
                 bg_prefixes: List[str], fg_seq_lengths: List[int],
                 bg_seq_lengths: List[int], max_extension: int = 70000,
                 uniformity_weight: float = 0.0,
                 reaction_temp: Optional[float] = None,
                 tm_weight: float = 0.0,
                 dimer_penalty: float = 0.0,
                 max_dimer_bp: int = 4,
                 conditions: Optional['ReactionConditions'] = None,
                 mechanistic_weight: float = 0.0,
                 template_gc: float = 0.5):
        """
        Initialize optimizer.

        Args:
            position_cache: PositionCache with primer positions
            fg_prefixes: Target genome prefixes
            bg_prefixes: Background genome prefixes
            fg_seq_lengths: Target genome lengths
            bg_seq_lengths: Background genome lengths
            max_extension: Max Phi29 extension (default 70 kb)
            uniformity_weight: Weight for coverage uniformity in scoring (0.0-1.0)
            reaction_temp: Reaction temperature for Tm-based scoring (None = disabled)
            tm_weight: Weight for Tm-based scoring (0.0-1.0)
            dimer_penalty: Penalty weight for primer-primer dimers (0.0-1.0)
            max_dimer_bp: Maximum complementary base pairs to consider as dimer
            conditions: ReactionConditions for mechanistic model (None = disabled)
            mechanistic_weight: Weight for mechanistic scoring (0.0-1.0)
            template_gc: Template genome GC content for accessibility calculation
        """
        self.cache = position_cache
        self.fg_prefixes = fg_prefixes
        self.bg_prefixes = bg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_seq_lengths = bg_seq_lengths
        self.max_extension = max_extension
        self.uniformity_weight = uniformity_weight
        self.reaction_temp = reaction_temp
        self.tm_weight = tm_weight
        self.dimer_penalty = dimer_penalty
        self.max_dimer_bp = max_dimer_bp
        # Cache for primer Tm values
        self._tm_cache: Dict[str, float] = {}

        # Mechanistic model integration
        self.conditions = conditions
        self.mechanistic_weight = mechanistic_weight
        self.template_gc = template_gc
        self.mech_model = None

        if conditions is not None and mechanistic_weight > 0:
            from neoswga.core.mechanistic_model import MechanisticModel
            self.mech_model = MechanisticModel(conditions)
            logger.info(f"Mechanistic model enabled (weight={mechanistic_weight:.2f})")

    def _get_primer_tm(self, primer: str) -> float:
        """Get effective primer Tm (cached).

        When :attr:`conditions` is set we use the canonical nearest-neighbor
        Tm plus additive correction from :meth:`ReactionConditions.calculate_tm_correction`
        so additives (DMSO, betaine, trehalose, formamide, ethanol, urea, TMAC)
        change Tm-weighted edge scoring. Without conditions we fall back to the
        legacy `melting_temp.temp()` estimate to preserve existing behaviour.
        """
        if primer in self._tm_cache:
            return self._tm_cache[primer]

        tm: float
        if self.conditions is not None:
            try:
                from neoswga.core.thermodynamics import calculate_tm_with_salt
                na = getattr(self.conditions, 'na_conc', 50.0)
                gc_count = sum(1 for b in primer if b in 'GC')
                primer_len = len(primer)
                gc_fraction = gc_count / primer_len if primer_len else 0.5
                tm = calculate_tm_with_salt(primer, na_conc=na)
                tm += self.conditions.calculate_tm_correction(
                    gc_content=gc_fraction, primer_length=primer_len,
                )
            except Exception:
                tm = calculate_primer_tm(primer)
        else:
            tm = calculate_primer_tm(primer)
        self._tm_cache[primer] = tm
        return tm

    def _calculate_tm_score(self, primer: str) -> float:
        """
        Calculate Tm-based scoring factor.

        Returns 0-1 where 1 = optimal Tm for reaction, 0 = poor Tm.
        Optimal Tm is ~5C above reaction temperature for primer binding.
        """
        if self.reaction_temp is None:
            return 1.0  # No Tm scoring if no reaction temp set

        tm = self._get_primer_tm(primer)
        optimal_tm = self.reaction_temp + 5  # Optimal annealing ~5C above reaction

        # Gaussian decay: efficiency drops as |Tm - optimal| increases
        delta = abs(tm - optimal_tm)
        return math.exp(-0.05 * delta**2)  # 50% at ~3.7C difference

    def _calculate_set_dimer_penalty(self, primer: str, current_set: List[str]) -> float:
        """
        Calculate dimer penalty for adding primer to current set.

        Returns 0-1 where 0 = no dimers, 1 = dimers with all existing primers.
        """
        if not current_set or self.dimer_penalty <= 0:
            return 0.0

        # Check self-dimer
        self_dimer = calculate_dimer_score(primer, primer, self.max_dimer_bp)

        # Check cross-dimers with existing primers
        cross_dimers = 0.0
        for existing in current_set:
            cross_dimers += calculate_dimer_score(primer, existing, self.max_dimer_bp)

        # Normalize by set size + 1 (for self)
        total_penalty = (self_dimer + cross_dimers) / (len(current_set) + 1)
        return total_penalty

    def _calculate_mechanistic_score(self, primer: str) -> float:
        """
        Calculate mechanistic model scoring factor.

        Returns 0-1 where 1 = optimal amplification conditions, 0 = poor conditions.
        Uses the four-pathway mechanistic model to predict amplification efficiency.
        """
        if self.mech_model is None:
            return 1.0

        effects = self.mech_model.calculate_effects(primer, self.template_gc)
        return effects.predicted_amplification_factor

    def optimize_greedy(self, candidates: List[str], num_primers: int = 10) -> List[str]:
        """
        Greedy network optimization.

        At each step, add primer that maximally increases target network
        connectivity while minimizing background network connectivity.
        Optionally weights coverage uniformity if uniformity_weight > 0.

        Args:
            candidates: Candidate primer list
            num_primers: Number of primers to select

        Returns:
            Selected primer set
        """
        logger.info(f"Optimizing primer set: {len(candidates)} candidates -> {num_primers} primers")
        if self.uniformity_weight > 0:
            logger.info(f"  Uniformity weight: {self.uniformity_weight:.2f}")
        if self.tm_weight > 0 and self.reaction_temp is not None:
            logger.info(f"  Tm weight: {self.tm_weight:.2f} (reaction temp: {self.reaction_temp}C)")
        if self.dimer_penalty > 0:
            logger.info(f"  Dimer penalty: {self.dimer_penalty:.2f} (max bp: {self.max_dimer_bp})")
        if self.mech_model is not None:
            logger.info(f"  Mechanistic weight: {self.mechanistic_weight:.2f} (template GC: {self.template_gc:.2f})")

        selected = []
        current_fg_network = AmplificationNetwork(self.max_extension)
        current_bg_network = AmplificationNetwork(self.max_extension)

        for iteration in range(num_primers):
            logger.info(f"Iteration {iteration+1}/{num_primers}")

            best_primer = None
            best_score = -float('inf')

            # Ensure spatial indices are clean before the candidate loop so
            # that _simulate_add_primer can save/restore them cheaply instead
            # of re-sorting the full binding_sites list for every candidate.
            current_fg_network._rebuild_spatial_index()
            current_bg_network._rebuild_spatial_index()

            # Try each remaining candidate
            for primer in candidates:
                if primer in selected:
                    continue

                # Evaluate adding this primer (base network score)
                score = self._evaluate_primer_addition(
                    primer, selected, current_fg_network, current_bg_network,
                    uniformity_weight=self.uniformity_weight
                )

                # Apply Tm weighting if enabled
                if self.tm_weight > 0 and self.reaction_temp is not None:
                    tm_factor = self._calculate_tm_score(primer)
                    # Blend: (1 - tm_weight) * base + tm_weight * tm_factor * base
                    score = score * (1.0 - self.tm_weight + self.tm_weight * tm_factor)

                # Apply dimer penalty if enabled
                if self.dimer_penalty > 0:
                    dimer_factor = self._calculate_set_dimer_penalty(primer, selected)
                    # Reduce score by dimer penalty
                    score = score * (1.0 - self.dimer_penalty * dimer_factor)

                # Apply mechanistic model weighting if enabled
                if self.mech_model is not None and self.mechanistic_weight > 0:
                    mech_factor = self._calculate_mechanistic_score(primer)
                    # Blend: boost score by mechanistic amplification factor
                    score = score * (1.0 - self.mechanistic_weight + self.mechanistic_weight * mech_factor)

                if score > best_score:
                    best_score = score
                    best_primer = primer

            if best_primer is None:
                logger.warning("No more primers can be added")
                break

            # Add best primer
            selected.append(best_primer)
            logger.info(f"  Selected: {best_primer} (score={best_score:.3f})")

            # Update networks
            self._update_network(current_fg_network, best_primer, self.fg_prefixes, 'target')
            self._update_network(current_bg_network, best_primer, self.bg_prefixes, 'background')

        # Final statistics
        fg_stats = current_fg_network.get_statistics()
        bg_stats = current_bg_network.get_statistics()

        logger.info(f"Final target network: {fg_stats['largest_component']} largest component, "
                   f"{fg_stats['predicted_amplification']:.0f}× predicted amplification")
        logger.info(f"Final background network: {bg_stats['avg_component_size']:.1f} avg component size")

        return selected

    def _evaluate_primer_addition(self, primer: str, current_set: List[str],
                                  fg_network: AmplificationNetwork,
                                  bg_network: AmplificationNetwork,
                                  uniformity_weight: float = 0.0) -> float:
        """
        Score of adding primer to current set.

        Base score = (target_network_improvement) / (background_network_penalty + 1)

        With uniformity weighting:
        Final score = (1 - uniformity_weight) * base_score + uniformity_weight * uniformity_improvement

        Args:
            primer: Primer to evaluate
            current_set: Current selected primers
            fg_network: Target genome network
            bg_network: Background genome network
            uniformity_weight: Weight for coverage uniformity (0.0-1.0)

        Returns:
            Combined score
        """
        # Simulate adding primer to target
        test_fg = self._simulate_add_primer(fg_network, primer, self.fg_prefixes)
        fg_improvement = test_fg.largest_component_size() - fg_network.largest_component_size()

        # Simulate adding primer to background
        test_bg = self._simulate_add_primer(bg_network, primer, self.bg_prefixes)
        bg_penalty = test_bg.average_component_size() - bg_network.average_component_size()

        # Base score: maximize target connectivity, minimize background connectivity
        base_score = fg_improvement / (bg_penalty + 1.0)

        # If uniformity weight is 0, return base score (backward compatible)
        if uniformity_weight <= 0.0:
            return base_score

        # Calculate uniformity improvement
        # Lower CV is better, so improvement = old_CV - new_CV
        genome_length = sum(self.fg_seq_lengths)
        old_cv = fg_network.coverage_uniformity(genome_length) if fg_network.binding_sites else float('inf')
        new_cv = test_fg.coverage_uniformity(genome_length)

        # Convert CV improvement to score (positive is better)
        # Normalize: CV typically 0-2, so improvement of 0.5 is significant
        if old_cv == float('inf'):
            uniformity_improvement = 1.0  # First primer always improves uniformity
        else:
            uniformity_improvement = max(0.0, min(1.0, (old_cv - new_cv) + 0.5))

        # Combine scores
        final_score = (1.0 - uniformity_weight) * base_score + uniformity_weight * uniformity_improvement

        return final_score

    def _simulate_add_primer(self, network: AmplificationNetwork,
                            primer: str, prefixes: List[str]) -> AmplificationNetwork:
        """
        Simulate adding primer to network using in-place add/remove.

        Adds primer sites and edges to the network, computes metrics,
        then removes them to restore the original state. Avoids the
        O(n+m) graph copy that dominated runtime in the previous approach.

        Returns a lightweight proxy with the computed metrics.
        """
        # Save original state for restoration (including spatial index so we
        # can restore it directly rather than re-sorting for every candidate)
        original_bs_len = len(network.binding_sites)
        saved_positions = network._sorted_positions
        saved_sites = network._sorted_sites
        saved_dirty = network._index_dirty

        # Add new primer sites in-place
        all_new_sites = []
        for prefix in prefixes:
            fw_positions = self.cache.get_positions(prefix, primer, 'forward')
            rv_positions = self.cache.get_positions(prefix, primer, 'reverse')

            new_sites_fw = network.add_primer_sites(primer, fw_positions, '+')
            new_sites_rv = network.add_primer_sites(primer, rv_positions, '-')
            all_new_sites.extend(new_sites_fw)
            all_new_sites.extend(new_sites_rv)

        # Build edges for new sites only
        network.add_edges_for_sites(all_new_sites)

        # Compute all metrics from a single connected_components traversal
        has_sites = bool(network.binding_sites)
        n_nodes = len(network.graph)
        if n_nodes == 0:
            largest = 0
            avg = 0.0
            uniformity = float('inf')
        else:
            components = list(nx.connected_components(network.graph))
            comp_sizes = [len(c) for c in components]
            largest = max(comp_sizes)
            avg = sum(comp_sizes) / len(comp_sizes)

            if len(components) < 2:
                uniformity = 0.0
            else:
                comp_positions = [np.mean([s.position for s in c]) for c in components]
                mean_pos = np.mean(comp_positions)
                uniformity = np.std(comp_positions) / mean_pos if mean_pos != 0 else 0.0

        # Restore original state: remove new nodes (and their edges)
        network.graph.remove_nodes_from(all_new_sites)
        network.binding_sites = network.binding_sites[:original_bs_len]
        # Restore the saved spatial index instead of marking dirty (avoids
        # O(n log n) re-sort for every candidate in the greedy loop)
        network._sorted_positions = saved_positions
        network._sorted_sites = saved_sites
        network._index_dirty = saved_dirty

        # Return lightweight proxy with cached metrics
        return _SimulatedNetwork(largest, avg, uniformity, has_sites)

    def _update_network(self, network: AmplificationNetwork, primer: str,
                       prefixes: List[str], label: str):
        """
        Update network with new primer using incremental edge building.

        Args:
            network: Network to update in place
            primer: Primer sequence to add
            prefixes: Genome prefixes to get positions from
            label: Label for logging ('target' or 'background')
        """
        all_new_sites = []
        for prefix in prefixes:
            fw_positions = self.cache.get_positions(prefix, primer, 'forward')
            rv_positions = self.cache.get_positions(prefix, primer, 'reverse')

            new_sites_fw = network.add_primer_sites(primer, fw_positions, '+')
            new_sites_rv = network.add_primer_sites(primer, rv_positions, '-')
            all_new_sites.extend(new_sites_fw)
            all_new_sites.extend(new_sites_rv)

        # Use incremental edge building
        network.add_edges_for_sites(all_new_sites)

    def score_primer_set(self, primers: List[str],
                         primer_tms: Optional[Dict[str, float]] = None,
                         reaction_temp: Optional[float] = None) -> Dict:
        """
        Comprehensive scoring of primer set with improved enrichment calculation.

        Args:
            primers: List of primer sequences
            primer_tms: Optional dict mapping primer to Tm for weighted scoring
            reaction_temp: Optional reaction temperature for weighted scoring

        Returns:
            dict with the following keys:

            - ``primers``: list of primer sequences evaluated.
            - ``num_primers``: number of primers in the set.
            - ``target_largest_component``: size of the largest connected
              component in the foreground amplification network.
            - ``target_connectivity``: algebraic connectivity of the
              foreground network.
            - ``target_amplification``: predicted foreground amplification
              fold (Tm-weighted when primer_tms is provided).
            - ``target_amplification_lower``: lower confidence bound on
              foreground amplification.
            - ``target_amplification_upper``: upper confidence bound on
              foreground amplification.
            - ``target_gini``: Gini coefficient of foreground binding
              positions (0 = uniform, 1 = concentrated).
            - ``background_avg_component``: mean connected-component size
              in the background network.
            - ``background_amplification``: predicted background
              amplification fold.
            - ``background_amplification_lower``: lower confidence bound
              on background amplification.
            - ``background_amplification_upper``: upper confidence bound
              on background amplification.
            - ``enrichment``: foreground / background amplification ratio.
            - ``enrichment_lower_ci``: conservative (worst-case) enrichment
              bound.
            - ``enrichment_upper_ci``: optimistic (best-case) enrichment
              bound.
            - ``score``: overall score (equal to ``enrichment``).

            When a mechanistic model is configured, these additional keys
            are included:

            - ``avg_processivity_factor``: mean polymerase processivity
              factor across primers.
            - ``avg_accessibility``: mean template accessibility factor.
            - ``avg_amplification_factor``: mean predicted mechanistic
              amplification factor.
            - ``min_amplification_factor``: minimum mechanistic
              amplification factor across primers.
            - ``avg_effective_binding``: mean effective primer binding
              rate.
            - ``avg_stability_factor``: mean enzyme stability factor.
        """
        # Build networks
        fg_network = self._build_network(primers, self.fg_prefixes)
        bg_network = self._build_network(primers, self.bg_prefixes)

        fg_stats = fg_network.get_statistics()
        bg_stats = bg_network.get_statistics()

        # Get amplification predictions
        fg_amplification = fg_stats['predicted_amplification']
        bg_amplification = bg_stats['predicted_amplification']

        # Use Tm-weighted predictions if available
        if primer_tms and reaction_temp:
            fg_amplification = fg_network.predict_amplification_fold_weighted(
                primer_tms, reaction_temp
            )
            bg_amplification = bg_network.predict_amplification_fold_weighted(
                primer_tms, reaction_temp
            )

        # Improved enrichment calculation:
        # - If background amplification is negligible (<1), enrichment equals
        #   foreground amplification (effectively infinite selectivity)
        # - Otherwise, compute ratio without arbitrary +1 offset
        if bg_amplification < 1.0:
            enrichment = fg_amplification  # No significant background
        else:
            enrichment = fg_amplification / bg_amplification

        # Get confidence intervals
        fg_lower = fg_stats.get('amplification_lower_ci', fg_amplification * 0.5)
        fg_upper = fg_stats.get('amplification_upper_ci', fg_amplification * 2.0)
        bg_lower = bg_stats.get('amplification_lower_ci', max(0.5, bg_amplification * 0.5))
        bg_upper = bg_stats.get('amplification_upper_ci', bg_amplification * 2.0)

        # Conservative enrichment bounds (worst case: low fg, high bg)
        enrichment_lower = fg_lower / max(bg_upper, 1.0)
        # Optimistic bounds (best case: high fg, low bg)
        enrichment_upper = fg_upper / max(bg_lower, 1.0)

        # Calculate mechanistic metrics if model is available
        mech_metrics = {}
        if self.mech_model is not None:
            effects_list = [
                self.mech_model.calculate_effects(p, self.template_gc)
                for p in primers
            ]
            if effects_list:
                import numpy as np
                mech_metrics = {
                    'avg_processivity_factor': float(np.mean([e.processivity_factor for e in effects_list])),
                    'avg_accessibility': float(np.mean([e.accessibility_factor for e in effects_list])),
                    'avg_amplification_factor': float(np.mean([e.predicted_amplification_factor for e in effects_list])),
                    'min_amplification_factor': float(min(e.predicted_amplification_factor for e in effects_list)),
                    'avg_effective_binding': float(np.mean([e.effective_binding_rate for e in effects_list])),
                    'avg_stability_factor': float(np.mean([e.stability_factor for e in effects_list])),
                }

        result = {
            'primers': primers,
            'num_primers': len(primers),
            'target_largest_component': fg_stats['largest_component'],
            'target_connectivity': fg_stats['connectivity'],
            'target_amplification': fg_amplification,
            'target_amplification_lower': fg_lower,
            'target_amplification_upper': fg_upper,
            'target_gini': fg_stats.get('gini_coefficient', 0),
            'background_avg_component': bg_stats['avg_component_size'],
            'background_amplification': bg_amplification,
            'background_amplification_lower': bg_lower,
            'background_amplification_upper': bg_upper,
            'enrichment': enrichment,
            'enrichment_lower_ci': enrichment_lower,
            'enrichment_upper_ci': enrichment_upper,
            'score': enrichment,  # Overall score
        }

        # Add mechanistic metrics if available
        result.update(mech_metrics)

        # Strand alternation analysis
        if self.cache is not None and hasattr(self.cache, 'compute_strand_alternation_stats'):
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                try:
                    strand_stats = self.cache.compute_strand_alternation_stats(
                        prefix, primers, length
                    )
                    result['strand_alternation_score'] = strand_stats['strand_alternation_score']
                    result['strand_coverage_ratio'] = strand_stats['strand_coverage_ratio']
                    result['strand_alternation_gap_mean'] = strand_stats['strand_alternation_gap_mean']
                    result['longest_same_strand_run'] = strand_stats['longest_same_strand_run']
                    break
                except Exception as e:
                    logger.debug(f"Strand analysis skipped: {e}")

        return result

    def _build_network(self, primers: List[str], prefixes: List[str]) -> AmplificationNetwork:
        """Build complete network for primer set"""
        network = AmplificationNetwork(self.max_extension)

        for primer in primers:
            for prefix in prefixes:
                fw_positions = self.cache.get_positions(prefix, primer, 'forward')
                rv_positions = self.cache.get_positions(prefix, primer, 'reverse')

                network.add_primer_sites(primer, fw_positions, '+')
                network.add_primer_sites(primer, rv_positions, '-')

        network.build_edges()
        return network


def benchmark_network_vs_ratio(position_cache, fg_prefixes, bg_prefixes,
                              fg_seq_lengths, bg_seq_lengths, candidates):
    """
    Compare network-based vs. ratio-based primer selection.

    Demonstrates superiority of network approach.
    """
    # Network-based selection
    network_opt = NetworkOptimizer(position_cache, fg_prefixes, bg_prefixes,
                                   fg_seq_lengths, bg_seq_lengths)
    network_primers = network_opt.optimize_greedy(candidates, num_primers=10)
    network_score = network_opt.score_primer_set(network_primers)

    # Ratio-based selection (current algorithm)
    ratio_primers = []
    for primer in candidates:
        fg_count = sum(len(position_cache.get_positions(p, primer)) for p in fg_prefixes)
        bg_count = sum(len(position_cache.get_positions(p, primer)) for p in bg_prefixes)
        ratio = fg_count / (bg_count + 1)
        ratio_primers.append((ratio, primer))

    ratio_primers.sort(reverse=True)
    ratio_primers = [p for _, p in ratio_primers[:10]]
    ratio_score = network_opt.score_primer_set(ratio_primers)

    logger.info("=== Network-based selection ===")
    logger.info(f"Primers: {network_primers}")
    logger.info(f"Target amplification: {network_score['target_amplification']:.0f}x")
    logger.info(f"Background amplification: {network_score['background_amplification']:.0f}x")
    logger.info(f"Enrichment: {network_score['enrichment']:.1f}x")

    logger.info("=== Ratio-based selection (current) ===")
    logger.info(f"Primers: {ratio_primers}")
    logger.info(f"Target amplification: {ratio_score['target_amplification']:.0f}x")
    logger.info(f"Background amplification: {ratio_score['background_amplification']:.0f}x")
    logger.info(f"Enrichment: {ratio_score['enrichment']:.1f}x")

    logger.info("=== Improvement ===")
    improvement = network_score['enrichment'] / ratio_score['enrichment']
    logger.info(f"Network method {improvement:.1f}x better enrichment")


# =============================================================================
# Factory Registration - BaseOptimizer Interface
# =============================================================================

from neoswga.core.base_optimizer import (
    BaseOptimizer, OptimizationResult, OptimizationStatus,
    PrimerSetMetrics, OptimizerConfig
)
from neoswga.core.optimizer_factory import OptimizerFactory


@OptimizerFactory.register('network', aliases=['network-optimizer', 'tm-weighted'])
class NetworkBaseOptimizer(BaseOptimizer):
    """
    Network optimizer implementing BaseOptimizer interface.

    Uses graph connectivity for Tm-weighted primer selection.
    """

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
        **kwargs
    ):
        super().__init__(
            position_cache, fg_prefixes, fg_seq_lengths,
            bg_prefixes, bg_seq_lengths, config,
            conditions=conditions,
        )
        self._network = NetworkOptimizer(
            position_cache=position_cache,
            fg_prefixes=fg_prefixes,
            bg_prefixes=bg_prefixes or [],
            fg_seq_lengths=fg_seq_lengths,
            bg_seq_lengths=bg_seq_lengths or [],
            max_extension=kwargs.get('max_extension', 70000),
            uniformity_weight=kwargs.get('uniformity_weight', 0.0),
            reaction_temp=kwargs.get('reaction_temp'),
            tm_weight=kwargs.get('tm_weight', 0.0),
            dimer_penalty=kwargs.get('dimer_penalty', 0.0),
            conditions=kwargs.get('conditions'),
            mechanistic_weight=kwargs.get('mechanistic_weight', 0.0),
            template_gc=kwargs.get('template_gc', 0.5),
        )

    @property
    def name(self) -> str:
        """Optimizer identifier for logging and factory registration."""
        return "network"

    @property
    def description(self) -> str:
        """One-line summary of the optimization strategy."""
        return "Network-based optimizer with Tm weighting and dimer penalty"

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        **kwargs
    ) -> OptimizationResult:
        """Run network-based optimization."""
        candidates = self._validate_candidates(candidates)
        target = target_size or self.config.target_set_size

        if self.config.verbose:
            logger.info(f"Running network optimization: {len(candidates)} candidates")

        try:
            primers = self._network.optimize_greedy(
                candidates=candidates,
                num_primers=target,
            )

            metrics = self.compute_metrics(primers)

            return OptimizationResult(
                primers=tuple(primers),
                score=metrics.fg_coverage,
                status=OptimizationStatus.SUCCESS if primers else OptimizationStatus.NO_CONVERGENCE,
                metrics=metrics,
                iterations=1,
                optimizer_name=self.name,
            )

        except Exception as e:
            logger.error(f"Network optimization failed: {e}")
            return OptimizationResult.failure(self.name, str(e))
