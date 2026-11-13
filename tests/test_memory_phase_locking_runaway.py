"""The search for a stable cycle abandons a start whose trajectory runs away (phase_locking.find_cycle, _solve). Found on
corpus cells such as the one-loop running of three gauge couplings, dg/dt = b g^3, which reaches a pole in finite time:
near the pole the LSODA step fell to about 1e-15 ("t + h = t"), solve_ivp kept every accepted step and the process grew
to tens of gigabytes; elsewhere a NaN in the analytic Jacobian reached LAPACK and crashed the process."""
from __future__ import annotations

import time
import tracemalloc

import numpy as np
import pytest

from fieldbridge.memory import phase_locking as pl
from fieldbridge.memory import spec as mspec

BASE = {"schema": "fieldbridge-memory/1", "kind": "equations", "field": "test", "question": "test", "assumptions": [],
        "provenance": {"source": "test", "contributor": "test"}, "noise": 0.0, "dt": 0.01}

RUNNING = {
    **BASE, "name": "one-loop running of three couplings (test)",
    "carrier": {"kind": "euclid", "name": "couplings", "variables": ["x", "y", "z"], "scale": 0.2},
    "parameters": {"n": 3.0},
    "drift": {"x": "(2*n - 9)*x**3", "y": "(2*n - 5)*y**3", "z": "(2*n + 3/5)*z**3"},
    "control": {"name": "n", "range": [2.9, 3.1]},
}

STIFF = {
    **BASE, "name": "stiff relaxation (test)",
    "carrier": {"kind": "euclid", "name": "plane", "variables": ["x", "y"], "scale": 1},
    "parameters": {"k": 1e4},
    "drift": {"x": "-k*(x - y)", "y": "-y"},
    "control": {"name": "k", "range": [1e3, 1e5]},
}

GEODESIC = {
    **BASE, "name": "Eguchi-Hanson geodesic flow (test, arXiv:2302.02620)",
    "carrier": {"kind": "euclid", "name": "radius, radial momentum, angular momenta",
                "variables": ["rho", "p_rho", "M1", "M2", "M3"], "scale": 1},
    "parameters": {"ga": 1.0},
    "drift": {"rho": "(rho**2-ga**2)*p_rho/rho",
              "p_rho": "(-p_rho**2*(rho**2+ga**2)/rho**2 + (M1**2+M2**2)/rho**2 + M3**2*(rho**2+ga**2)/(rho**2-ga**2)**2)/2",
              "M1": "-ga**2*M2*M3/(rho*(rho**2-ga**2))", "M2": "ga**2*M3*M1/(rho*(rho**2-ga**2))", "M3": "0"},
    "control": {"name": "ga", "range": [0.1, 2.0]},
}

VDP = {
    **BASE, "name": "van der Pol (test)",
    "carrier": {"kind": "euclid", "name": "plane", "variables": ["x", "v"], "scale": 1},
    "parameters": {"mu": 1.0},
    "drift": {"x": "v", "v": "mu*(1 - x**2)*v - x"},
    "control": {"name": "mu", "range": [0.5, 2.0]},
}


def test_a_start_that_reaches_a_pole_in_finite_time_is_abandoned_with_bounded_memory():
    # y and z grow as (t* - t)^(-1/2) from every start; without the budget one start ran for minutes and its chunk held
    # hundreds of megabytes of steps
    real = mspec.load(RUNNING)
    tracemalloc.start()
    t0 = time.perf_counter()
    try:
        cyc = pl.find_cycle(real, {}, np.random.default_rng(0))
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert cyc is None
    assert time.perf_counter() - t0 < 60.0
    assert peak < 200 * 2 ** 20
    rng = np.random.default_rng(0)
    with pytest.raises(pl._Abandoned, match="step collapsed"):
        pl._solve(real, {}, (0.0, 100.0), real.carrier.sample(rng, 1)[0], rtol=1e-9, atol=2e-12, budget=True)


def test_a_non_finite_jacobian_ends_the_run_before_lapack_factors_it(monkeypatch):
    real = mspec.load(STIFF)
    sol = pl._solve(real, {}, (0.0, 10.0), [1.0, 1.0], rtol=1e-9, atol=1e-11)
    assert sol.status == 0 and np.allclose(sol.y[:, -1], np.exp(-10.0), rtol=1e-3)
    monkeypatch.setattr(pl._Compiled, "jac", lambda self, q, pv: np.full((len(q), len(q)), np.nan))
    sol = pl._solve(real, {}, (0.0, 10.0), [1.0, 1.0], rtol=1e-9, atol=1e-11)
    assert sol.status < 0 or not np.all(np.isfinite(sol.y))
    try:
        sol = pl._solve(real, {}, (0.0, 10.0), [1.0, 1.0], rtol=1e-9, atol=1e-11, budget=True)
    except pl._Abandoned as e:
        assert "Jacobian" in str(e)
    else:
        assert sol.status < 0


def test_a_chunk_over_its_step_budget_abandons_the_start(monkeypatch):
    real = mspec.load(VDP)
    assert pl.find_cycle(real, {}, np.random.default_rng(0)) is not None
    monkeypatch.setattr(pl, "CHUNK_STEPS", 50)  # a chunk of 100 time units of the cycle takes thousands of steps
    assert pl.find_cycle(real, {}, np.random.default_rng(0)) is None


def test_a_collapsed_step_that_cannot_grow_again_still_ends_the_run():
    # from this start rho reaches 0, where the drift has a 1/rho singularity, at t = 0.491; the watch finds the collapsed
    # step and sets the drift to zero, but LSODA's step stays near 2e-18, below the rounding of t, so t never moves and
    # solve_ivp stored every return: over a million drift calls at one t, gigabytes within minutes on the cluster
    real = mspec.load(GEODESIC)
    q = real.carrier.sample(np.random.default_rng(6), 1)[0]
    tracemalloc.start()
    t0 = time.perf_counter()
    try:
        with np.errstate(all="ignore"), pytest.raises(pl._Abandoned, match="step collapsed"):
            pl._solve(real, {}, (0.0, 100.0), q, rtol=1e-9, atol=1e-11, budget=True)
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert time.perf_counter() - t0 < 30.0
    assert peak < 200 * 2 ** 20
