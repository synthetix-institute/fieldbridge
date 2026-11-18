"""The derivatives of the normal form in the control where the drift is not defined on one side (correction 28 of the
discovery engine).

normal_form takes the derivatives of the drift and of its Jacobian in the control as central differences over v -+ dv,
and the eigenvalues of the Jacobian at v + dv for the slope of the critical rate. Where the control is an exponent and a
variable is negative, the drift at a non-integer exponent is not defined: the Jacobian at v + dv was not finite, and the
eigenvalues raised "Array must not contain infs or NaNs" (a corpus model whose exponent the scan took to m = 70). The
derivatives are now taken on the side where the drift is defined, and where it is defined on neither they are zero."""
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


# the exponent as the control: at m = 2 the drift is defined at x = -0.5, at m = 2 -+ dv it is not
EXPONENT = _spec("power with the exponent as control", ["x", "y"], {"m": 2.0}, {"x": "y", "y": "-x**m - y"}, "m",
                 [1.5, 2.5])
# defined only below p = 1: at p = 1 - dv / 2 the drift at p + dv is not defined, at p - dv it is
ONE_SIDED = _spec("square root of the control", ["x"], {"p": 0.0}, {"x": "sqrt(1 - p) - x"}, "p", [0.0, 1.0])
# the supercritical pitchfork, with every value finite
PITCHFORK = _spec("pitchfork", ["x"], {"r": 0.0}, {"x": "r*x - x**3"}, "r", [-0.5, 0.5])


def test_a_drift_undefined_on_both_sides_of_the_control_gives_a_normal_form():
    real = spec.load(EXPONENT)
    q = np.array([-0.5, 0.0])
    assert np.all(np.isfinite(an.jacobian(real, q, {"m": 2.0})))
    assert not np.all(np.isfinite(an.jacobian(real, q, {"m": 2.0 + 1e-5})))
    out = an.normal_form(real, q, "m", 2.0)
    assert isinstance(out["kind"], str)
    assert out["kappa_slope"] == 0.0 or np.isfinite(out["kappa_slope"])


def test_a_drift_undefined_on_one_side_is_differenced_on_the_other():
    real = spec.load(ONE_SIDED)
    dv = 1e-5 * 1.0
    v = 1.0 - 0.5 * dv
    q = np.array([np.sqrt(1.0 - v)])
    out = an.normal_form(real, q, "p", v)
    # the drift along x has the bias d/dp sqrt(1 - p) = -1 / (2 sqrt(1 - p)), taken backwards from p - dv
    exact = -1.0 / (2.0 * np.sqrt(1.0 - v))
    backward = (np.sqrt(1.0 - v) - np.sqrt(1.0 - v + dv)) / dv
    assert np.isfinite(out["bias_along_mode"])
    assert abs(abs(out["bias_along_mode"]) - abs(backward)) < 1e-6 * abs(backward)
    assert abs(backward - exact) < 0.5 * abs(exact)


def test_where_every_value_is_finite_the_normal_form_is_unchanged():
    out = an.normal_form(spec.load(PITCHFORK), np.array([0.0]), "r", 0.0)
    assert out["kind"].startswith("supercritical pitchfork")
    assert abs(out["a3"] + 1.0) < 1e-6
    assert abs(out["kappa_slope"] + 1.0) < 1e-6
    assert out["bias"] < 1e-9
