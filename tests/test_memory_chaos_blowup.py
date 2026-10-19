"""Chaotic loss of the initial state and blow-up in finite time (fieldbridge.memory.chaos_blowup).

The Lorenz system is chaotic with lambda_1 close to 0.906 (Sprott, 2003), and the time for a separation to grow obeys
the predictability law. Blow-up: x' = x^2 + 1 grows as (t* - t)^-1; x' = y, y' = x^2 (x'' = x^2) has degree 2 but
grows as (t* - t)^-3 in its fastest component, which only the dominant balance predicts; a damped oscillator does not
escape, and the logistic-free linear growth x' = x escapes without blowing up.
"""
import numpy as np
import pytest

pytest.importorskip("sympy")
pytest.importorskip("scipy")

from fieldbridge.memory import chaos_blowup as cb
from fieldbridge.memory import spec as S


def model(variables, drift, params=None, scale=1.0):
    return S.load({"schema": "fieldbridge-memory/1", "kind": "equations", "name": "test", "question": "test",
                   "assumptions": [], "provenance": "test", "noise": 0.0,
                   "carrier": {"kind": "euclid", "variables": variables, "scale": scale},
                   "parameters": params or {}, "drift": drift})


LORENZ = model(["x", "y", "z"], {"x": "s*(y - x)", "y": "x*(r - z) - y", "z": "x*y - b*z"},
               {"s": 10.0, "r": 28.0, "b": 8.0 / 3}, scale=10.0)


def test_the_lorenz_system_loses_its_initial_state_at_its_lyapunov_exponent():
    r = cb.chaos(LORENZ, np.random.default_rng(0), starts=3, transient=100.0, T=300.0, span=200.0)
    assert r["chaotic"]
    assert abs(r["lyapunov"]["median"] - 0.906) < 0.05
    assert r["law_holds"] and abs(r["predictability"]["law"] - 1.0) < 0.2


def test_a_damped_oscillator_is_not_chaotic():
    r = cb.chaos(model(["x", "v"], {"x": "v", "v": "-x - 0.5*v"}), np.random.default_rng(0), starts=3, transient=50.0, T=60.0, span=40.0)
    assert not r["chaotic"]


def test_the_balance_of_the_dominant_terms_predicts_the_growth_where_the_degree_does_not():
    quad = cb.blowup(model(["x"], {"x": "x**2 + 1"}), np.random.default_rng(1), starts=3)
    assert abs(quad["growth_exponent"] - 1.0) < 0.05 and quad["predicted_by_balance"]["growth"] == pytest.approx(1.0)
    assert quad["escape_time_scaling"]["weighted"] and quad["escape_time_scaling"]["slope"] == pytest.approx(-1.0, abs=0.05)
    p1 = cb.blowup(model(["x", "y"], {"x": "y", "y": "x**2"}), np.random.default_rng(1), starts=3)
    assert p1["degree"] == 2
    assert abs(p1["growth_exponent"] - 3.0) < 0.3                       # not 1/(2 - 1) = 1
    assert p1["predicted_by_balance"]["p"] == pytest.approx([2.0, 3.0])  # x ~ tau^-2, y ~ tau^-3
    assert p1["predicted_by_balance"]["c"] == pytest.approx([6.0, 12.0], rel=1e-4)
    # starts scaled with the weights (2, 3) escape after a time proportional to 1/s; the subdominant term still
    # shifts the slope by a few per cent at these amplitudes
    assert p1["escape_time_scaling"]["weighted"] and p1["escape_time_scaling"]["slope"] == pytest.approx(-1.0, abs=0.1)


def test_growth_without_blow_up_and_no_escape_are_told_apart():
    lin = cb.blowup(model(["x"], {"x": "x"}), np.random.default_rng(2), starts=3)
    assert lin["verdict"].startswith("growth without blow-up")
    damped = cb.blowup(model(["x", "v"], {"x": "v", "v": "-x - 0.5*v"}), np.random.default_rng(2), starts=3)
    assert damped["verdict"] == "no start escapes"
