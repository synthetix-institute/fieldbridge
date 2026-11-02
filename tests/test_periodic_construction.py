"""Tests for the pulse-parity toggle and amount-and-phase constructions."""
import json
import copy
from pathlib import Path

import numpy as np
import pytest

from fieldbridge.memory import periodic_construction as pc

TOGGLE_SPEC = {
    "schema": "fieldbridge-memory/1",
    "question": "pulse-parity toggle test",
    "assumptions": ["test"],
    "provenance": {"source": "test", "kind": "authored example"},
    "kind": "equations",
    "name": "bistable_adiabatic_closure",
    "carrier": {"kind": "euclid", "name": "x", "variables": ["x"]},
    "parameters": {"mu": 1.5},
    "drift": {"x": "mu*x - x**3 - 0.45*x"},
    "control": {"name": "mu", "range": [-1.0, 2.0]},
    "material": [],
    "noise": 0.015,
    "adiabatic_closure": {"variable": "z", "kappa": 0.45, "unreduced_drift": "mu*x - x**3"},
}

NETWORK_SPEC = {
    "schema": "fieldbridge-memory/1",
    "question": "cascaded network test",
    "assumptions": ["test"],
    "provenance": {"source": "arXiv:0903.5248", "kind": "extracted"},
    "kind": "equations",
    "name": "model_0903_5248",
    "carrier": {"kind": "orthant", "name": "c", "variables": ["x_1", "x_2", "x_3", "x_4"]},
    "parameters": {"c": 1.0, "b": 3.3, "c_": 2.0},
    "drift": {
        "x_1": "-x_1*(b - c_*x_1*x_2 + 1) + 1",
        "x_2": "x_1*(b - c_*x_1*x_2)",
        "x_3": "-x_3*(-c*x_3*x_4 + x_2 + 1) + 1",
        "x_4": "x_3*(-c*x_3*x_4 + x_2)",
    },
    "control": {"name": "c_", "range": [1.0, 3.0]},
    "material": ["c", "b"],
    "noise": 0.01,
}


def test_one_dimensional_obstruction():
    out = pc.check_1d_toggle_obstruction(TOGGLE_SPEC)
    assert out["obstructed"] and out["odd_drift"] and out["zero_invariant"]


def test_closure_promotion_selects_slow_negative_feedback():
    res = pc.derive_pulse_parity_toggle(TOGGLE_SPEC)
    assert res["status"] == "verified_construction"
    status = {r["id"]: r["status"] for r in res["closure_promotions"]}
    assert status == {
        "multiplicative_coupling": "rejected_symbolically",
        "positive_additive_feedback": "rejected_symbolically",
        "negative_additive_fast_closure": "rejected_dynamically",
        "negative_additive_slow_closure": "selected",
    }
    sel = res["selected_construction"]
    assert max(sel["period_two_multiplier_moduli"]) == pytest.approx(0.015548, abs=2e-4)
    q0, q1 = np.array(sel["period_two_states"])
    assert np.allclose(q0, -q1, atol=1e-8)
    assert abs(q0[0]) == pytest.approx(1.087145, abs=1e-4)
    assert not res["stochastic_verification"]["fast_closure_control"]["certified_alternation"]


def test_cascade_hopf_thresholds():
    out = pc.analyze_orthant_hopf_structure(NETWORK_SPEC)
    assert out["upstream_hopf_control"] == pytest.approx(3.3 - 1.0)
    assert out["downstream_hopf_control"] == pytest.approx(3.3 / 2.0)


def test_phase_candidates_preserve_circle():
    base, cands = pc.generate_amount_phase_candidates(NETWORK_SPEC)
    assert base["parameters"]["g_add"] == 0.0
    assert len(cands) == 4 and all(c["unit_circle_invariant"] for c in cands)
    assert "g_add*(x_1*c_ph)" in cands[0]["spec"]["drift"]["x_1"]


def test_line_obstruction_is_not_applied_to_a_circle():
    spec = copy.deepcopy(TOGGLE_SPEC)
    spec["carrier"]["kind"] = "torus"
    assert pc.check_1d_toggle_obstruction(spec)["obstructed"] is None
    assert pc.derive_pulse_parity_toggle(spec)["status"] == "unsupported_source"


def test_closure_must_reproduce_source():
    spec = copy.deepcopy(TOGGLE_SPEC)
    spec["drift"]["x"] = "mu*x-7*x**3-0.45*x"
    assert pc.derive_pulse_parity_toggle(spec)["status"] == "unsupported_source"


def test_changed_cubic_is_used_by_the_compiled_equations():
    spec = copy.deepcopy(TOGGLE_SPEC)
    spec["drift"]["x"] = "mu*x-7*x**3-0.45*x"
    spec["adiabatic_closure"]["unreduced_drift"] = "mu*x-7*x**3"
    candidate = pc.enumerate_closure_promotions(spec)[-1]
    compiled = pc._compile_promotion(spec, candidate)
    assert compiled["real"].F(np.array([1., .45]))[0] == pytest.approx(1.5-7-.45)
    jac = compiled["jacobian"](np.array([1., .45]), 1.5)
    assert jac[0, 0] == pytest.approx(1.5-21)


def test_promotion_is_invariant_under_variable_renaming():
    spec = copy.deepcopy(TOGGLE_SPEC)
    spec["carrier"]["variables"] = ["amp"]
    spec["drift"] = {"amp": "mu*amp-amp**3-0.45*amp"}
    spec["adiabatic_closure"] = {"variable": "feedback", "kappa": .45,
                                  "unreduced_drift": "mu*amp-amp**3"}
    compiled = pc._compile_promotion(spec, pc.enumerate_closure_promotions(spec)[-1])
    assert compiled["real"].F(np.array([1., .45]))[0] == pytest.approx(.05)


def test_changed_source_changes_the_calculated_orbit(monkeypatch):
    seen = []
    def stochastic(*args, **kwargs):
        seen.append(kwargs["compiled"]["spec"])
        return {"certified_alternation": args[4] > 1}
    monkeypatch.setattr(pc, "_stochastic_toggle_check", stochastic)
    spec = copy.deepcopy(TOGGLE_SPEC)
    spec["drift"]["x"] = "mu*x-7*x**3-0.45*x"
    spec["adiabatic_closure"]["unreduced_drift"] = "mu*x-7*x**3"
    result = pc.derive_pulse_parity_toggle(spec)
    assert result["status"] == "verified_construction"
    selected = result["selected_construction"]
    assert abs(selected["period_two_states"][0][0]) == pytest.approx(1.0871454829/np.sqrt(7), abs=1e-7)
    assert selected["period_two_residual"] < 1e-8
    assert all("7*x**3" in s["drift"]["x"] for s in seen)
