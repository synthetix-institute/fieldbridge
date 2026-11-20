"""The budget of explicit steps of a relaxation (correction 31 of the discovery engine).

A relaxation stayed explicit up to 1e7 steps of the fastest rate. At 1e-4 s a step of 48 starts that is about 1000 s,
more than the engine's time limit, and the card relaxes more than once: a corpus epidemic of ten variables needed 6.3e6
steps and reached the time limit in its first relaxation. A relaxation above RELAX_BUDGET = 3e6 steps now goes to the
implicit method after the explicit prefix, as one above 1e7 did. Below 3e6 nothing changes."""
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


def _two_rates(K):
    # a variable slaved at the rate K to a slow one that relaxes at the rate 0.5: the fastest rate over the control range
    # is 2 K, so the whole time of 300 needs about 2000 K explicit steps; the slow variable settles near t = 20, after the
    # prefix of 1e5 steps when K = 2.5e3 (t = 6)
    return _spec("slow relaxation with a fast variable", ["x", "y"], {"K": K},
                 {"x": "-0.5*x + 0.01*y", "y": "-K*(y - x)"}, "K", [0.5 * K, 2 * K])


def _relax(real, monkeypatch):
    called = []
    implicit = an._implicit_runs

    def spy(*a, **k):
        called.append(True)
        return implicit(*a, **k)

    monkeypatch.setattr(an, "_implicit_runs", spy)
    q, speed = an.relax(real, np.array([[1.0, 1.0], [-0.8, 0.2]]))
    return q, speed, bool(called)


def test_a_relaxation_of_5e6_steps_goes_implicit(monkeypatch):
    real = spec.load(_two_rates(2.5e3))
    assert 3e6 < 300.0 / an.safe_dt(real) < 1e7
    q, speed, implicit = _relax(real, monkeypatch)
    assert implicit
    assert np.all(speed < 1e-4) and np.all(np.abs(q) < 1e-2)


def test_a_relaxation_of_2e6_steps_stays_explicit(monkeypatch):
    real = spec.load(_two_rates(1e3))
    assert 300.0 / an.safe_dt(real) < 3e6
    q, speed, implicit = _relax(real, monkeypatch)
    assert not implicit
    assert np.all(speed < 1e-4)
