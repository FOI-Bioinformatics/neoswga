"""`--min-fg-bg-ratio` has to reach the background pre-filter.

Its help text says "Minimum foreground/background binding site ratio. Higher
values = more selective". It was read in exactly two places, both of them
optional side paths -- `recommend_set_size` under `--auto-size` and
`select_from_frontier` under `--show-frontier` -- and never forwarded into the
optimizer. So `neoswga optimize --min-fg-bg-ratio 20` on its own applied no
selectivity constraint at all, silently, and the flag was not in
`_common.UNIMPLEMENTED_OPTIONS` either, so nothing warned.

The knob it should drive already existed and had no route to it:
`run_optimization` reads `kwargs["bg_min_ratio"]` as the pre-filter threshold
and `run_step4` never passed it.
"""

import pytest


class _Args:
    """The attributes `_step4_optimizer_kwargs` reads."""

    use_position_cache = True
    use_background_filter = False
    optimization_method = "hybrid"
    quiet = True


def test_the_flag_is_forwarded_as_the_prefilter_threshold():
    from neoswga.cli.pipeline import _step4_optimizer_kwargs

    args = _Args()
    args.min_fg_bg_ratio = 20.0

    assert _step4_optimizer_kwargs(args)["bg_min_ratio"] == 20.0


def test_an_unset_flag_forwards_nothing_to_override_the_default():
    from neoswga.cli.pipeline import _step4_optimizer_kwargs

    assert _step4_optimizer_kwargs(_Args())["bg_min_ratio"] is None


@pytest.mark.parametrize("forwarded,expected", [(None, 1.0), (20.0, 20.0)])
def test_run_optimization_hands_the_value_to_the_prefilter(monkeypatch, forwarded, expected):
    """Driven through the real `run_optimization`, not a copy of its line.

    None must resolve to the default rather than being passed through:
    `kwargs.get("bg_min_ratio", 1.0)` returns None when the key is present and
    None, which is exactly what the CLI now sends for an unset flag, so the
    default has to come from `or` rather than from `.get`.
    """
    from neoswga.core import unified_optimizer as uo

    seen = {}

    class _Stop(Exception):
        """Ends the run once the pre-filter has been reached; everything past
        it needs real position data and is not what this test is about."""

    def spy(cache, candidates, fg_prefixes, bg_prefixes, min_ratio, max_removal_fraction, verbose):
        seen["min_ratio"] = min_ratio
        seen["max_removal"] = max_removal_fraction
        raise _Stop

    monkeypatch.setattr(uo, "_prefilter_by_background", spy)

    class _Cache:
        def __init__(self, *a, **k):
            pass

    monkeypatch.setattr(uo, "PositionCache", _Cache)

    with pytest.raises(_Stop):
        uo.run_optimization(
            method="hybrid",
            candidates=["ACCACAGATAGC", "GTTGTAGATGGA"],
            fg_prefixes=["fg"],
            fg_seq_lengths=[10_000],
            bg_prefixes=["bg"],
            bg_seq_lengths=[10_000],
            target_size=2,
            verbose=False,
            bg_min_ratio=forwarded,
        )

    assert seen["min_ratio"] == expected
    assert seen["max_removal"] == 0.20
