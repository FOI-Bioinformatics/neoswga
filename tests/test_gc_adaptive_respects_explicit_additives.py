"""GC-adaptive defaults do not overwrite an additive the user set explicitly.

Finding E2. The guards at pipeline.py:559 and :566 tested `== 0.0`, which
cannot distinguish an explicit zero from an unset default, so a params.json
saying `"betaine_m": 0.0` was silently run at 1.0 M. The k-mer branch above them
already tests key presence in `parameter._json_data`; these now do the same.
"""

import neoswga.core.pipeline as pipeline_mod
from neoswga.core import parameter


class _AdaptiveParams:
    """The subset of GCAdaptiveStrategy output the additive branches read."""

    betaine_concentration = 1.0
    dmso_concentration = 5.0


def test_explicit_zero_betaine_is_preserved(monkeypatch):
    monkeypatch.setattr(parameter, "_json_data", {"betaine_m": 0.0}, raising=False)
    monkeypatch.setattr(parameter, "betaine_m", 0.0, raising=False)
    pipeline_mod._apply_adaptive_additives(_AdaptiveParams())
    assert parameter.betaine_m == 0.0


def test_explicit_zero_dmso_is_preserved(monkeypatch):
    monkeypatch.setattr(parameter, "_json_data", {"dmso_percent": 0.0}, raising=False)
    monkeypatch.setattr(parameter, "dmso_percent", 0.0, raising=False)
    pipeline_mod._apply_adaptive_additives(_AdaptiveParams())
    assert parameter.dmso_percent == 0.0


def test_unset_betaine_still_takes_the_recommendation(monkeypatch):
    monkeypatch.setattr(parameter, "_json_data", {}, raising=False)
    monkeypatch.setattr(parameter, "betaine_m", 0.0, raising=False)
    pipeline_mod._apply_adaptive_additives(_AdaptiveParams())
    assert parameter.betaine_m == 1.0


def test_unset_dmso_still_takes_the_recommendation(monkeypatch):
    monkeypatch.setattr(parameter, "_json_data", {}, raising=False)
    monkeypatch.setattr(parameter, "dmso_percent", 0.0, raising=False)
    pipeline_mod._apply_adaptive_additives(_AdaptiveParams())
    assert parameter.dmso_percent == 5.0
