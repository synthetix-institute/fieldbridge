"""The grid of a noiseless trajectory whose count of steps exceeds a 64-bit integer (correction 33 of the discovery engine).

A noiseless trajectory longer than the explicit budget is integrated by the implicit method at the times of a grid
built from the count of steps. Where the step of the fastest rate is of order 1e-36, as when the control range of a card
reaches an exponent of 600, a trajectory of 200 time units has 1e38 steps: the counts multiplied as integers gave an array
of Python integers, the grid an array of objects, and the implicit runs raised "ufunc 'isfinite' not supported for the
input types". Beyond a 64-bit integer the counts are now multiplied in floating point. Below it nothing changes."""
from __future__ import annotations

import numpy as np

from fieldbridge.memory import analysis as an
from fieldbridge.memory import spec


def _spec(name, variables, parameters, drift, control, control_range, scale=1.0):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": name, "question": "what does it store?",
            "assumptions": ["test model"], "provenance": {"source": "test"},
            "carrier": {"kind": "euclid", "name": "state", "variables": variables, "scale": scale},
            "parameters": parameters, "drift": drift, "control": {"name": control, "range": control_range},
            "noise": 0.01, "dt": 0.01}


ROTATION = _spec("slow rotation", ["x", "y"], {"w": 1.0}, {"x": "-w*y", "y": "w*x"}, "w", [0.5, 1.5])


def _trajectory(step):
    real = spec.load(ROTATION)
    real._safe_dt = step  # the step that safe_dt would give for a fastest rate of 0.3 / step
    return an.integrate(real, np.array([[1.0, 0.0]]), 200.0, 0.0, np.random.default_rng(0),
                        record=lambda x: x[0].copy(), n_record=40)


def test_a_count_of_steps_beyond_64_bits_gives_a_trajectory():
    q, t, v = _trajectory(2e-36)
    assert t.dtype == float and len(t) == 40
    assert np.all(np.isfinite(q))
    assert abs(t[-1] - 200.0) < 1e-6
    assert abs(np.linalg.norm(q) - 1.0) < 1e-2  # the rotation keeps its radius


def test_a_count_within_64_bits_gives_the_same_grid():
    step = 1e-6
    q, t, v = _trajectory(step)
    steps = int(round(200.0 / step))
    every = steps // 40
    assert np.array_equal(t, np.arange(1, steps // every + 1) * every * step)
