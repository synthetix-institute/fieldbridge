"""The reading rule of the two-mode surface (examples/ptsb): a share is a product ratio only where the later well is
bound. Statistical thresholds only."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples" / "ptsb"))
from two_modes import TwoModes, fate, well_depths  # noqa: E402


def test_the_z_well_gets_shallow_with_the_delay():
    assert well_depths(TwoModes(c=1.5, D=0.1))["z"] == pytest.approx(0.054, abs=0.002)
    assert well_depths(TwoModes(c=1.5, D=0.2))["z"] == pytest.approx(0.011, abs=0.002)


def test_a_share_in_a_shallow_well_is_a_transient():
    r = fate(TwoModes(c=1.5, D=0.2), 1.0, 0.05, n=1500, times=(40,), seed=7)
    assert not r["z_bound"]
    assert r["z_share_first_entry"] > 0.15 and r["z_share"][40] < 0.06
    bound = fate(TwoModes(c=1.5, D=0.05), 1.0, 0.01, n=1500, times=(20,), seed=7)
    assert bound["z_bound"]
    assert bound["z_share"][20] == pytest.approx(bound["z_share_first_entry"], abs=0.05)
