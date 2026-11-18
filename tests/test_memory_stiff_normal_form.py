"""The kind of a write point with a fast variable (correction 27 of the discovery engine).

normal_form reads the reduced drift as vanishing, a degenerate write point, when its quadratic and cubic terms over the
window of the fit are below 1e-6 of the rates times the window. With the largest rate of the spectrum as that rate, a
variable slaved at the rate 1e6 put the bound at the order of the reduced drift, and every fold and pitchfork carried
by a fast variable was read as degenerate. The rate is now the largest the declared step resolves, at most 0.3 / dt, as
in safe_dt. A drift that does vanish along the critical direction, as on a line of states, is degenerate as before."""
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


# a transcritical point whose quadratic term is carried by a variable slaved at the rate 1e6: reduced x' = r x - x^2
TRANSCRITICAL = _spec("transcritical point with a fast mediator", ["x", "y"], {"r": 0.0, "K": 1e6},
                      {"x": "r*x - x*y", "y": "-K*(y - x)"}, "r", [-0.5, 0.5])
# a pitchfork whose cubic term is the fast variable itself: reduced x' = r x - x^3
PITCHFORK = _spec("pitchfork with a fast cubic", ["x", "y"], {"r": 0.0, "K": 1e6},
                  {"x": "r*x - y", "y": "-K*(y - x**3)"}, "r", [-0.5, 0.5])
# a fast exchange between two compartments with a slow loss mu: at mu = 0 every equal pair is a state
EXCHANGE = _spec("fast exchange with a slow loss", ["a", "b"], {"mu": 0.0, "K": 1e6},
                 {"a": "-K*(a - b) - mu*a", "b": "K*(a - b) - mu*b"}, "mu", [-0.5, 0.5])
# the reduced pitchfork itself, with no fast rate
REDUCED = _spec("pitchfork", ["x"], {"r": 0.0}, {"x": "r*x - x**3"}, "r", [-0.5, 0.5])


def _kind(s, q, param="r", v=0.0):
    return an.normal_form(spec.load(s), np.array(q, float), param, v)["kind"]


def test_a_transcritical_point_with_a_fast_mediator_is_transcritical():
    assert _kind(TRANSCRITICAL, [0.0, 0.0]).startswith("transcritical")


def test_a_pitchfork_with_a_fast_cubic_is_a_supercritical_pitchfork():
    assert _kind(PITCHFORK, [0.0, 0.0]).startswith("supercritical pitchfork")


def test_a_line_of_states_with_a_fast_exchange_is_still_degenerate():
    assert _kind(EXCHANGE, [0.5, 0.5], "mu", 0.0).startswith("degenerate")


def test_without_a_fast_rate_the_kind_is_unchanged():
    assert _kind(REDUCED, [0.0]).startswith("supercritical pitchfork")
    real = spec.load(REDUCED)
    out = an.normal_form(real, np.array([0.0]), "r", 0.0)
    scale = max(abs(float(np.linalg.eigvals(an.jacobian(real, np.array([0.0]), {"r": 0.0}))[0])), 1e-12)
    assert scale <= 0.3 / real.dt  # the bound is the one of the code before
    assert abs(out["a3"] + 1.0) < 1e-6
