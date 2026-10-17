"""The swept-write law across the noise ratio Gamma = D|a3|/r: pi^{1/4} for weak noise, an activated law for strong
noise, the same function for continuous and first-order switches."""
import json
from math import pi
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from fieldbridge.memory import noise_ratio as nr
from fieldbridge.memory import spec

ROOT = Path(__file__).resolve().parents[1]
K = pi ** 0.25


def test_the_law_starts_at_kondepudi_and_falls_monotonically():
    assert nr.write_constant(1e-3) == pytest.approx(K) and nr.write_constant(1.0) == pytest.approx(K)
    grid = np.logspace(0, 2.5, 30)
    for first_order in (False, True):
        c = [nr.write_constant(g, first_order) for g in grid]
        assert all(b <= a + 1e-12 for a, b in zip(c, c[1:]))
    for g in nr.GAMMAS:
        assert abs(nr.write_constant(g, True) / nr.write_constant(g) - 1) < 0.04
    assert [nr.regime(g) for g in (0.05, 1.0, 10.0, 30.0)] == [nr.LINEAR, nr.LINEAR, nr.CROSSOVER, nr.ACTIVATED]


def test_the_activated_limit_is_an_upper_estimate_at_large_noise_ratio():
    for g in (30.0, 100.0):
        act = nr.activated_constant(g)
        assert 1.0 < act / nr.write_constant(g) < 1.15


@pytest.mark.parametrize("gamma,first_order", [(0.1, False), (10.0, False), (10.0, True)])
def test_the_normal_form_reproduces_the_table(gamma, first_order):
    r = nr.normal_form_constant(gamma, 6000, np.random.default_rng(7), first_order)
    assert abs(r["constant"] - nr.write_constant(gamma, first_order)) < 4 * r["constant_stderr"] + 0.02


def test_a_first_order_switch_decides_before_its_bifurcation_on_the_normal_form():
    r = nr.normal_form_sweeps(3.0, 0.0, 4000, np.random.default_rng(3), first_order=True)
    assert r["decided_before"] > 0.95
    assert nr.normal_form_sweeps(3.0, 0.0, 2000, np.random.default_rng(3))["decided_before"] == 0.0


def test_a_continuous_switch_follows_kondepudi_at_small_noise_ratio():
    real = spec.load(ROOT / "examples/memory/pitchfork.json")
    nf = nr.symmetric_point(real)
    rep = nr.sweeps(real, nf, np.random.default_rng(5), gammas=(0.05,), n_traj=600)
    row = rep["rows"][0]
    assert rep["first_order"] is False and row["regime"] == nr.LINEAR
    assert abs(row["constant"] - K) < 4 * row["constant_stderr"]


def test_a_first_order_switch_has_the_same_constant_and_decides_early():
    real = spec.load(ROOT / "examples/noise_ratio/first_order_landau.json")
    nf = nr.symmetric_point(real)
    assert nf is not None and nf["a3"] > 0
    rep = nr.sweeps(real, nf, np.random.default_rng(5), gammas=(1.0,), n_traj=400)
    row = rep["rows"][0]
    assert rep["first_order"] is True
    assert abs(row["constant"] - nr.write_constant(1.0, True)) < 4 * row["constant_stderr"]
    assert row["decided_before"] > 0.6


def test_the_command_writes_the_law(tmp_path):
    from fieldbridge.memory.cli import cmd_noise_ratio
    args = SimpleNamespace(spec=None, law=True, out_dir=str(tmp_path), seed=1, gammas="1", q=5e-4, trajectories=10)
    assert cmd_noise_ratio(args) == 0
    report = json.loads((tmp_path / "noise_ratio.json").read_text())
    assert report["law"]["continuous"][0] == pytest.approx(K)
    assert "activated" in (tmp_path / "noise_ratio.md").read_text()
