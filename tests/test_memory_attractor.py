"""The kind of motion of a realization without a stable state (memory.attractor): chaos, a cycle, a quasi-periodic
motion, a drift, convergence, escape. Statistical and numerical verdicts on known systems; no exponent is pinned."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from fieldbridge.memory import attractor as at
from fieldbridge.memory import spec

ROOT = Path(__file__).resolve().parents[1]


def _spec(name, drift, parameters, variables, carrier="euclid", scale=1.0, dt=0.01, noise=0.01):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": name, "field": "test",
            "question": "what does the motion settle on?", "assumptions": ["test model"], "provenance": {"source": "test"},
            "carrier": {"kind": carrier, "name": "state", "variables": variables, "scale": scale},
            "parameters": parameters, "drift": drift, "noise": noise, "dt": dt}


def test_the_lorenz_system_is_chaotic_and_the_van_der_pol_oscillator_a_cycle():
    lorenz = spec.load(ROOT / "examples" / "dynamics" / "lorenz.json")
    out = at.attractor(lorenz, np.random.default_rng(0), steps_cap=100_000, spec_steps=200_000)
    assert out["kind"] == "chaotic" and out["lyapunov"][0] > out["tolerance"], out
    vdp = spec.load(ROOT / "examples" / "memory" / "oscillators" / "van_der_pol.json")
    out = at.attractor(vdp, np.random.default_rng(0), steps_cap=100_000, spec_steps=200_000)
    assert out["kind"] == "cycle", out
    assert abs(out["lyapunov"][0]) <= out["tolerance"] and out["lyapunov"][1] < -out["tolerance"]
    assert out["spectrum"]["independent"] == 1


def test_two_uncoupled_oscillators_with_incommensurate_frequencies_are_quasi_periodic():
    s = _spec("two rotations", {"x1": "x1*(1 - x1**2 - y1**2) - 1.0*y1", "y1": "y1*(1 - x1**2 - y1**2) + 1.0*x1",
                                "x2": "x2*(1 - x2**2 - y2**2) - 1.4142135*y2", "y2": "y2*(1 - x2**2 - y2**2) + 1.4142135*x2"},
              {"w": 1.0}, ["x1", "y1", "x2", "y2"])
    out = at.attractor(spec.load(s), np.random.default_rng(0), steps_cap=100_000, spec_steps=200_000)
    assert out["kind"] == "quasi-periodic", out
    assert out["spectrum"]["independent"] == 2
    assert abs(out["spectrum"]["ratio"] - 2 ** 0.5) < 0.05 or abs(out["spectrum"]["ratio"] - 2 ** -0.5) < 0.05


def test_convergence_drift_and_escape_are_told_apart():
    damped = spec.load(_spec("damped", {"x": "-x", "y": "-2*y"}, {"k": 1.0}, ["x", "y"]))
    assert at.attractor(damped, np.random.default_rng(0), steps_cap=100_000)["kind"] == "converges"
    drift = spec.load(_spec("drift", {"x": "1 + 0*x"}, {"k": 1.0}, ["x"]))
    out = at.attractor(drift, np.random.default_rng(0), steps_cap=100_000)
    assert out["kind"] == "drift", out
    assert out["recurrence"]["crossings"] < at.REC_MIN
    blowup = spec.load(_spec("blow-up", {"x": "x**2 + 0.1"}, {"k": 1.0}, ["x"]))   # every start escapes within a few time units
    assert at.attractor(blowup, np.random.default_rng(0), steps_cap=100_000)["kind"] == "unbounded"
    # hyperbolic growth from a start near zero blows up only after thousands of time units: within the run it reads as
    # a drift or an escape, and the blow-up law (chaos_blowup) decides; the record runs that law for a drift
    growth = spec.load(ROOT / "examples" / "dynamics" / "hyperbolic_growth.json")
    assert at.attractor(growth, np.random.default_rng(0), steps_cap=50_000)["kind"] in ("drift", "unbounded")


def test_recurrence_and_frequencies_on_synthetic_records():
    t = np.arange(0, 400, 0.01)
    circle = np.stack([np.cos(t), np.sin(t)], axis=1)
    rec = at.recurrence(circle)
    assert rec["kind"] == "cycle" and rec["crossings"] >= at.REC_MIN
    torus = np.stack([np.cos(t) + 0.5 * np.cos(2 ** 0.5 * t), np.sin(t), np.sin(2 ** 0.5 * t)], axis=1)
    f = at.frequencies(torus, 0.01, np.random.default_rng(0))
    assert f["independent"] == 2
    line = np.stack([t, 0 * t], axis=1)
    assert at.recurrence(line)["kind"] == "drift"
