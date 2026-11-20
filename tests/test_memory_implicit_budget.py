"""A budget of calls of the drift for one implicit run of a relaxation (correction 30 of the discovery engine).

The relaxation integrates a start by backward differentiation when its explicit steps would exceed the budget
(correction 21). A run that creeps along a singular line of the drift, or follows a fast motion that never settles, took
small but finite steps without end: the inspiral of a corpus model made 1.5e5 calls of the drift in 20 s and reached the
time limit of the card. A run of a relaxation now ends, as not settled, after IMPLICIT_CALLS calls of the drift. The stiff
relaxations of the development data need at most 484."""
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


# a fast rotation that never settles: over 0.1 time units BDF follows about 300 periods with some 7e4 calls of the drift
ROTATION = _spec("fast rotation", ["x", "y"], {"w": 2e4}, {"x": "-w*y", "y": "w*x"}, "w", [1e4, 3e4])
# a stiff relaxation onto a stable state: a variable slaved at the rate 1e6
STIFF = _spec("stiff relaxation", ["x", "y"], {"K": 1e6}, {"x": "-x + y", "y": "-K*(y - 0.5*x)"}, "K", [5e5, 2e6])


def test_a_run_that_does_not_settle_ends_after_its_budget():
    real = spec.load(ROTATION)
    counted = [0]
    drift = real.F

    def F(q, **kw):
        counted[0] += 1
        return drift(q, **kw)

    real.F = F
    q, _ = an._implicit_runs(real, np.array([[1.0, 0.0]]), {}, 0.1, tol=1e-5)
    assert not np.all(np.isfinite(q))   # not settled
    assert counted[0] < 21_000          # the budget, and the few calls of the drift made by the Jacobians


def test_a_stiff_relaxation_settles_within_the_budget():
    real = spec.load(STIFF)
    q, _ = an._implicit_runs(real, np.array([[1.0, 0.0], [-0.5, 0.3]]), {}, 300.0, tol=1e-5)
    assert np.all(np.isfinite(q))
    assert np.all(np.linalg.norm(q, axis=1) < 1e-3)
