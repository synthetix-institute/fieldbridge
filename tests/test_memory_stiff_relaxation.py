"""Relaxation of a stiff drift above the explicit step budget.

The relaxation behind stored_states integrates the drift by explicit Euler for 300 time units at the step of the
fastest rate, 0.3 / lambda_max. With rates of 1e5 and more against slow ones (a fast exchange, a large population, rate
constants of chemistry) it needs 1e8 to 1e17 steps and reached the engine's time limit. Above EXPLICIT_BUDGET steps each
start is now integrated by backward differentiation over the same time, with the same speed criterion, and so is a
noiseless trajectory; below the budget nothing changes."""
from __future__ import annotations

import time

import numpy as np
import pytest

from fieldbridge.memory import analysis as an
from fieldbridge.memory import spec


def _spec(name, variables, parameters, drift, control, control_range, carrier="euclid"):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": name, "question": "what does it store?",
            "assumptions": ["test model"], "provenance": {"source": "test"},
            "carrier": {"kind": carrier, "name": "state", "variables": variables, "scale": 1.0},
            "parameters": parameters, "drift": drift, "control": {"name": control, "range": control_range},
            "noise": 0.01, "dt": 0.01}


STIFF = _spec("fast exchange with a slow source", ["x", "y"], {"K": 1e7, "s": 1.0},
              {"x": "-K*(x - y)", "y": "s - y"}, "s", [0.5, 1.5])
MILD = _spec("slow exchange with a slow source", ["x", "y"], {"K": 10.0, "s": 1.0},
             {"x": "-K*(x - y)", "y": "s - y"}, "s", [0.5, 1.5])
QUICK = _spec("fast relaxation alone", ["x"], {"K": 1e7, "s": 1.0}, {"x": "-K*(x - s)"}, "s", [0.5, 1.5])


def test_a_stiff_relaxation_finds_its_one_state_without_the_explicit_steps():
    real = spec.load(STIFF)
    assert 300.0 / an.safe_dt(real) > an.EXPLICIT_BUDGET
    t0 = time.perf_counter()
    states, counts, spectra, unconverged, continuum = an.stored_states(real, np.random.default_rng(0), 48)
    assert time.perf_counter() - t0 < 120
    assert len(states) == 1 and unconverged == 0 and not continuum
    assert states[0] == pytest.approx([1.0, 1.0], abs=1e-4)


def test_below_the_budget_the_explicit_relaxation_is_unchanged(monkeypatch):
    real = spec.load(MILD)
    q0 = real.carrier.sample(np.random.default_rng(1), 8)
    q, speed = an.relax(real, q0)
    monkeypatch.setattr(an, "_implicit_runs", lambda *a, **k: (_ for _ in ()).throw(AssertionError("implicit")))
    q2, speed2 = an.relax(real, q0)
    assert np.array_equal(q, q2) and np.array_equal(speed, speed2)


def test_a_stiff_relaxation_that_settles_within_the_prefix_is_unchanged(monkeypatch):
    # x' = -1e7 (x - 1) settles in about a hundred Euler steps, although 300 time units would take 1e10 of them
    real = spec.load(QUICK)
    assert 300.0 / an.safe_dt(real) > an.EXPLICIT_BUDGET
    q0 = real.carrier.sample(np.random.default_rng(1), 8)
    monkeypatch.setattr(an, "_implicit_runs", lambda *a, **k: (_ for _ in ()).throw(AssertionError("implicit")))
    q, speed = an.relax(real, q0)
    assert np.all(speed < 1e-5) and q == pytest.approx(np.ones((8, 1)), abs=1e-9)


def test_a_noiseless_stiff_trajectory_is_sampled_as_the_explicit_one_would_be():
    real = spec.load(STIFF)
    q0 = np.array([[0.0, 0.0]])
    q, t, traj = an.integrate(real, q0, 20.0, 0.0, np.random.default_rng(0), record=lambda x: x[0].copy(), n_record=40)
    assert len(t) == 40 and t[-1] == pytest.approx(20.0, rel=1e-6)
    # the slow variable follows 1 - exp(-t), and the fast one follows it
    assert traj[:, 1] == pytest.approx(1.0 - np.exp(-t), abs=1e-4)
    assert traj[:, 0] == pytest.approx(traj[:, 1], abs=1e-4)
    assert q[0] == pytest.approx([1.0 - np.exp(-20.0)] * 2, abs=1e-4)
