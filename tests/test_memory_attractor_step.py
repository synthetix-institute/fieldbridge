"""The step of the attractor's Runge-Kutta run (attractor.attractor; correction 35 of the discovery engine).

The run takes fixed counts of fourth-order Runge-Kutta steps. Correction 34 bounded safe_dt by |a| / |l|^2 over the
stable eigenvalues, which explicit Euler needs at a focus that turns much faster than it decays. A run at that step
covers a time shorter by the same factor: in the engine, Lorenz-84 ran for 415 time units instead of 4000, and a corpus
model for 0.116 instead of 1663, whose torus then read a drift. Runge-Kutta keeps such a focus stable at the step of the
fastest rate, rate_dt, which the run takes."""
from __future__ import annotations

import numpy as np
import pytest

from fieldbridge.memory import analysis as an
from fieldbridge.memory import attractor as at
from fieldbridge.memory import spec

# the van der Pol oscillator driving a fast mode that the motion excites and that decays slowly: u'' = -W^2 u - c u' + x,
# whose rates -c / 2 +- i W bound the Euler step at about c / (2 W^2) = 1.25e-5, against 0.3 / W = 7.5e-3 for the
# fastest rate
VDP_FAST_MODE = {"schema": "fieldbridge-memory/1", "kind": "equations", "name": "van der Pol with a weakly damped fast mode",
                 "question": "what does it settle on?", "assumptions": ["test model"], "provenance": {"source": "test"},
                 "carrier": {"kind": "euclid", "name": "state", "variables": ["x", "y", "u", "v"], "scale": 2.0},
                 "parameters": {"mu": 1.0, "W": 40.0, "c": 0.04},
                 "drift": {"x": "y", "y": "mu*(1 - x**2)*y - x", "u": "v", "v": "-W**2*u - c*v + x"},
                 "control": {"name": "c", "range": [0.03, 0.05]}, "noise": 0.01, "dt": 0.01}


def test_the_two_steps_differ_by_the_euler_bound():
    real = spec.load(VDP_FAST_MODE)
    assert an.rate_dt(real) / an.safe_dt(real) > 100


def test_the_attractor_runs_at_the_runge_kutta_step():
    real = spec.load(VDP_FAST_MODE)
    cap = 3000
    out = at.attractor(real, np.random.default_rng(1), steps_cap=cap, spec_steps=20_000)
    assert out["capped"]
    # the run covers the cap times the step of the fastest rate; at the Euler step it covered about 0.04
    assert out["time"] == pytest.approx(cap * an.rate_dt(real), rel=1e-9)
    assert out["time"] > 20.0
