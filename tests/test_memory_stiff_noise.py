"""Noisy trajectories of a stiff drift beyond the explicit step budget (correction 25 of the discovery engine).

The hold time and the write trajectories integrate the drift with noise by Euler-Maruyama at the step of the fastest
rate. With a fast rate of 1e5 or more that is 1e8 steps and beyond, and the card reached the engine's time limit there.
Above EXPLICIT_BUDGET steps such a trajectory now takes the drift-implicit step (I - dt J)^-1 (dt f + noise) at the
declared step of the specification: a rate far above 1 / dt relaxes within a step, so the fast variable is slaved to
the slow ones, while the slow motion and its escape over a barrier are integrated as before. Below the budget nothing
changes."""
from __future__ import annotations

import time

import numpy as np
import pytest

from fieldbridge.memory import analysis as an
from fieldbridge.memory import spec


def _spec(name, variables, parameters, drift, noise=0.08):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": name, "question": "what does it store?",
            "assumptions": ["test model"], "provenance": {"source": "test"},
            "carrier": {"kind": "euclid", "name": "state", "variables": variables, "scale": 1.5},
            "parameters": parameters, "drift": drift, "noise": noise, "dt": 0.01}


# a double well in x with a variable y slaved to x at the rate K: with y = x the drift of x is x - x^3
STIFF = _spec("double well with a fast slaved variable", ["x", "y"], {"K": 1e7, "c": 1.0},
              {"x": "x - x**3 + c*(y - x)", "y": "-K*(y - x)"})
REDUCED = _spec("double well", ["x"], {}, {"x": "x - x**3"})
MILD = _spec("double well with a slow slaved variable", ["x", "y"], {"K": 3.0, "c": 1.0},
             {"x": "x - x**3 + c*(y - x)", "y": "-K*(y - x)"})


def _hold(s, n_traj=400):
    real = spec.load(s)
    states = [np.full(real.carrier.dim, v) for v in (-1.0, 1.0)]
    return an.hold_time(real, 0, states, real.noise, np.random.default_rng(3), n_traj=n_traj, t_max=300.0)


def _left(h, t):
    return float(np.asarray(h["fraction"])[np.searchsorted(h["t"], t)])


def test_the_stiff_escape_is_that_of_the_reduced_well():
    real = spec.load(STIFF)
    assert 300.0 / an.safe_dt(real) > an.EXPLICIT_BUDGET
    t0 = time.perf_counter()
    stiff = _hold(STIFF)
    assert time.perf_counter() - t0 < 120
    reduced = _hold(REDUCED)
    # Kramers: barrier 1/4 at D = 0.08, a mean escape time of order 1e2. The fractions of 400 trajectories still in the
    # first well differ by sampling only (standard deviation 0.025 each); the first time below 0.9 is too noisy to compare
    for t in (60.0, 120.0):
        assert _left(stiff, t) == pytest.approx(_left(reduced, t), abs=0.1)
    assert 0.3 < _left(stiff, 120.0) < 0.85


def test_a_noisy_trajectory_below_the_budget_is_unchanged():
    real = spec.load(MILD)
    assert 300.0 / an.safe_dt(real) < an.EXPLICIT_BUDGET
    q0 = np.full((5, 2), -1.0)
    q, t, v = an.integrate(real, q0.copy(), 20.0, real.noise, np.random.default_rng(7), record=lambda x: x[:, 0].mean())
    # the explicit Euler-Maruyama of the code before the correction
    rng, dt = np.random.default_rng(7), an.safe_dt(real)
    p, amp = q0.copy(), np.sqrt(2.0 * real.noise * dt)
    for _ in range(int(round(20.0 / dt))):
        p = real.carrier.wrap(p + dt * real.F(p) + amp * rng.standard_normal(p.shape))
    assert np.array_equal(q, p)


def test_a_fast_variable_is_slaved_within_one_step():
    # one drift-implicit step from y far from x: y lands on x to within 1 / (K dt)
    real = spec.load(STIFF)
    q = np.array([[0.5, -0.5]])
    kick = 0.01 * real.F(q)
    step = an._drift_implicit(real, q, {}, 0.01, kick)
    assert abs((q + step)[0, 1] - (q + step)[0, 0]) < 1e-4
