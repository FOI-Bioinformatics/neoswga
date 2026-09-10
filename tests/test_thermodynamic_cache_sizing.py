"""The two thermodynamic caches are sized from one named, measured constant.

The size was two bare literals, and a filter step reported holding 1,639 of the
million entries it was provisioned for. `lru_cache` allocates lazily, so the
figure is an eviction ceiling rather than a memory cost, and the value is right
for the largest shipped configuration -- but nothing in the code said either
thing, and a reader had no way to tell a measured number from an arbitrary one.
"""

import logging

import pytest

from neoswga.core import thermodynamics as thermo


def test_both_caches_are_sized_from_the_same_constant():
    assert (
        thermo.calculate_enthalpy_entropy_cached.cache_info().maxsize == thermo.THERMO_CACHE_MAXSIZE
    )
    assert (
        thermo.compute_free_energy_for_two_strings_cached.cache_info().maxsize
        == thermo.THERMO_CACHE_MAXSIZE
    )


def test_the_ceiling_covers_the_largest_shipped_candidate_pool():
    """369,459 candidates reach the scan in tests/validation/genomes.

    Measured on the plasmid example, step 2 makes about 1.9 enthalpy calls per
    candidate (22,175 calls for 11,803 candidates). Calls bound the number of
    distinct keys rather than counting them, since repeated calls on the same
    sequence collapse, so two keys per candidate is a ceiling on the demand and
    not a measurement of it. Sizing below that bound could introduce eviction on
    the runs that are already slowest.
    """
    assert thermo.THERMO_CACHE_MAXSIZE >= 2 * 369_459


def test_the_reported_capacity_is_the_configured_one(caplog):
    thermo.clear_thermodynamic_caches()
    thermo.calculate_enthalpy_entropy("ACGTACGTACGT")

    with caplog.at_level(logging.INFO, logger="neoswga.core.thermodynamics"):
        thermo.log_cache_stats("test")

    assert (
        f"/{thermo.THERMO_CACHE_MAXSIZE:,} entries" in caplog.text
    ), f"log_cache_stats did not report the configured ceiling; got {caplog.text!r}"
