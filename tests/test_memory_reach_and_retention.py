"""What a completed write law needs besides its moments: a sweep that comes to the write point, a reference that stays a state the sweep can
follow, a forcing counted once, a sign that survives saturation, and written states that keep the choice.

The numbers come from the equations: exact Gaussian moments with Jacobians written by hand, computed outside this package. The law's moments carry no
sampling error, and their gate is the second-order step of their integration."""
from __future__ import annotations

from math import pi, sqrt

import numpy as np
from fieldbridge.memory import analysis as an, completion as cm, construct, spec

RATE = 0.02
PI4 = pi ** 0.25


def real_of(variables, parameters, drift, control, lo=-1.0, hi=1.0, noise=0.001, carrier="euclid"):
    return spec.load({"schema": "fieldbridge-memory/1", "kind": "equations", "name": "test model", "field": "test", "question": "what the write follows", "assumptions": ["a test model"],
                      "provenance": {"source": "constructed for this test", "kind": "test"}, "carrier": {"kind": carrier, "name": "state", "variables": list(variables), "scale": 1.0},
                      "parameters": parameters, "drift": drift, "noise": noise, "control": {"name": control, "range": [lo, hi]}, "material": [k for k in parameters if k != control],
                      "observable": "the sign of the critical coordinate", "closure": "none"})


def law_at_zero(real, q=None, certify=False, samples=4000, traj=400, seed=0, **kw):
    q = np.zeros(real.carrier.dim) if q is None else np.asarray(q, float)
    nf = an.normal_form(real, q, real.control, 0.0)
    return cm.law_at(real, nf, np.random.default_rng(seed), RATE, None, (1.0,), samples, traj, certify, 0.05, **kw)


def moments(law):
    z = [r["mean_critical"][0] / sqrt(r["cov_critical"][0][0]) for r in law["rows"]]
    return PI4 * 0.5 * abs(z[0] - z[1]), 0.5 * (z[0] + z[1])


SHEAR = {"x": "p*(x - alpha*p) - (x - alpha*p)**3"}


def populations(mu):
    """x = u + s, y = u - s: u' = u (mu + u - u^2) keeps an empty state, and s crosses at p = 0 on the active state u0."""
    u, s = "(0.5*(x + y))", "(0.5*(x - y))"
    fu = f"{u}*({mu} + {u} - {u}**2)"
    fs = f"(p - g*(1 - {u}/u0))*{s} - 20*{s}**3"
    return {"x": f"{fu} + {fs}", "y": f"{fu} - ({fs})"}


def test_a_sweep_that_starts_on_another_state_has_no_law():
    # the active state ends in a fold at about p = -0.26, so at the start of the sweep only the empty state exists, and it stays stable
    mu0 = -0.15
    u0 = 0.5 * (1 + sqrt(1 + 4 * mu0))
    real = real_of(["x", "y"], {"p": 0.0, "g": 1.5, "u0": u0}, populations("(-0.1 + 0.25*(p - sqrt(p**2 + 0.04)))"), "p", carrier="orthant")
    law = law_at_zero(real, [u0, u0])
    assert law["status"].startswith("the sweep does not reach the write point")
    assert max(law["start"]) < 1e-3
    assert law["distance_at_crossing"] > 0.5


def test_a_sweep_that_lags_behind_the_branch_still_reaches_the_write_point():
    mu0 = -0.1
    u0 = 0.5 * (1 + sqrt(1 + 4 * mu0))
    real = real_of(["x", "y"], {"p": 0.0, "g": 1.5, "u0": u0}, populations("(-0.1 + 0.05*p)"), "p", carrier="orthant")
    law = law_at_zero(real, [u0, u0])
    assert law["status"] == "derived"
    assert law["lag_at_crossing"] > 1e-4
    # the rate is the slope of the critical rate along the branch, and nothing else changes it at first order
    assert abs(moments(law)[0] - 1.330) < 0.01


def test_a_flow_held_along_the_crossing_direction_is_not_followed_far_from_the_crossing():
    # the state is stable at every p < 0 through the coupling of the two modes, while with s held on its branch the relaxing mode has a rate
    # -g - e p, positive at the start of the sweep: holding it from there lets the reference run away
    real = real_of(["s", "w"], {"p": 0.0, "g": 0.3, "e": 1.2, "al": 0.3},
                   {"s": "p*s + p*(w - al*p) - 0.2*(w - al*p)**2*s - s**3", "w": "-g*(w - al*p) - p*s - e*p*(w - al*p)"}, "p")
    law = law_at_zero(real)
    assert law["status"] == "derived"
    assert law["lag_at_crossing"] < 0.5


def test_a_state_held_at_the_boundary_of_an_orthant_is_not_reflected():
    # the crossing direction R = (1, -1.3) and the relaxing one (1, -1) are nearly parallel, and the sweep starts on the boundary x = 0: a
    # reflection of the held state would multiply its excursion by (1.3 + 1) / (1.3 - 1) = 7.7 at every step
    bp, bq, y0 = 1.0, 1.3, 2.5
    w = f"(({bq}*x + (y - {y0}))/({bq} - {bp}))"
    s = f"(({bp}*x + (y - {y0}))/({bp} - {bq}))"
    fw, fs = f"-0.5*{w}", f"p*{s} - {s}**3"
    real = real_of(["x", "y"], {"p": 0.0}, {"x": f"{fw} + ({fs})", "y": f"-{bp}*({fw}) - {bq}*({fs})"}, "p", carrier="orthant")
    q = np.array([0.0, y0])
    nf = an.normal_form(real, q, "p", 0.0)
    R, L = cm.critical_space(real, nf, 1)
    sweep = cm.Sweep(real, nf, RATE)
    path = cm.reference_path(real, sweep, q, R, L, (R, L), [q + 10 * R[:, 0]], [1e-3], 0.01, 1e-3)
    assert np.all(np.isfinite(np.array(path["q"])))
    assert max(float(np.linalg.norm(x - q)) for x in path["q"]) < 1e-3
    assert path["ended_by"] == "growth"


def test_the_moments_follow_the_exact_ones_to_the_second_order_of_the_step():
    real = real_of(["s"], {"p": 0.0}, {"s": "p*s - s**3"}, "p")
    assert abs(moments(law_at_zero(real))[0] - 1.3303) < 5e-4          # the exact moments; an Euler step gives 1.3352


def test_the_forcing_of_the_sweep_is_counted_once():
    real = real_of(["x"], {"p": 0.0, "alpha": 0.3}, SHEAR, "p")
    law = law_at_zero(real, samples=64000)
    assert abs(abs(moments(law)[1]) - 0.6712) < 0.002                  # the exact offset, from the moments
    assert abs(abs(law["one_direction"]["law"][0]["offset"]) - 0.6712) < 0.03    # and from the probabilities; counted twice it is 0.728


def test_in_one_direction_a_saturated_sample_goes_by_its_sign():
    # no cubic term: the rotation saturates the growth, and the noise makes the samples at the mapping time larger than the written amplitude
    real = real_of(["x", "y"], {"p": 0.0, "g": 1.0, "K": 0.05}, {"x": "p*x - K*(x**2 + y**2)*y", "y": "-g*y + K*(x**2 + y**2)*x"}, "p", noise=0.05)
    law = law_at_zero(real, samples=16000)
    assert law["rows"][0]["regime"].startswith("saturated")
    c, _ = moments(law)
    assert abs(c - 1.3303) < 1e-3
    assert abs(law["one_direction"]["law"][0]["response_constant"] - c) < 0.05     # a flow from those samples gives about 0.78


def test_a_law_that_agrees_with_the_certificate_is_an_equilibrium_write_if_the_hold_exchanges_the_states(monkeypatch):
    real = real_of(["s"], {"p": 0.0}, {"s": "p*s - s**3"}, "p")
    monkeypatch.setattr(cm, "hold_check", lambda *a, **k: {"hold_time": 1.0, "switching_fraction": 0.5, "switching_by_row": [0.5, 0.5], "after_hold": []})
    law = law_at_zero(real, certify=True)
    assert law["certified"] is True
    assert law["write_class"] == "equilibrium write"
    assert "agrees" in law["obstruction_in_nonlinear_stage"]
    monkeypatch.setattr(cm, "hold_check", lambda *a, **k: {"hold_time": 1.0, "switching_fraction": 0.0, "switching_by_row": [0.0, 0.0], "after_hold": []})
    assert law_at_zero(real, certify=True)["write_class"] == "linear-stage write"


def test_a_failure_at_one_second_parameter_keeps_the_later_ones(monkeypatch):
    # a write point that needs a second parameter, unfolded along two of them; the first fails as the reference of version 2 did
    real = real_of(["s"], {"p": 0.0, "b": 1.0, "c": 1.0}, {"s": "p*s - s**3"}, "p")
    monkeypatch.setattr(cm.an, "locate_writes", lambda *a, **k: [{"q": np.zeros(1), "v": 0.0}])
    monkeypatch.setattr(cm.an, "normal_form", lambda *a, **k: {"kind": "saddle-node (fold)", "state": [0.0], "param": "p", "value": 0.0})
    monkeypatch.setattr(cm.construct, "obstruction", lambda nf: {"asymmetric": True, "terms": {}})

    def cancel(real_, nf, second):
        if second == "b":
            raise np.linalg.LinAlgError("Array must not contain infs or NaNs")
        return {"success": False}

    monkeypatch.setattr(cm.construct, "cancel_asymmetry", cancel)
    entry = cm.complete(real, np.random.default_rng(0), certify=False)["write_points"][0]
    assert entry["completion"] == "second parameter"
    assert [u["second"] for u in entry["unfoldings"]] == ["b", "c"]
    assert entry["unfoldings"][0]["law"]["status"].startswith("not derived")
    assert "law" not in entry


def test_the_followed_branch_stays_unwritten_through_the_crossing():
    # the unwritten state x = 0.5 T moves along the critical direction while the fictive temperature lags: at r = 0.01 Newton from the last point slid
    # onto a written state at the crossing, and the critical rate never turned positive (0.207 against 1.347)
    real = real_of(["x", "f"], {"T": 0.0, "B": 1.0, "G0": 0.6, "Ea": 3.0, "al": 0.5},
                   {"x": "-B*f*(x - al*T) - (x - al*T)**3", "f": "G0*exp(Ea*T)*(T - f)"}, "T", -0.8, 0.8)
    nf = an.normal_form(real, np.zeros(2), "T", 0.0)
    sweep = cm.Sweep(real, nf, 0.01)
    R, L = cm.critical_space(real, nf, 1)
    q0, how = cm.start_of_sweep(real, sweep, np.zeros(2), (R, L))
    path = cm.reference_path(real, sweep, q0, R, L, (R, L), [1e3 * R[:, 0]], [1e-3], 0.01, 1e-3)
    assert path["ended_by"] == "growth"
    t, B = np.array(path["t"]), np.array(path["branch"])
    after = t > sweep.t_cross + 5.0
    assert np.all(np.abs(B[after, 0] - 0.5 * np.array([sweep.value(x) for x in t[after]])) < 1e-3)    # on x = al T, not on a written state


def test_the_sweep_starts_on_the_unwritten_branch_where_the_write_point_lies_outside_its_basin():
    # w relaxes to 0.5 p only from within 0.27 of it; at the start the write point w = 0 lies 0.28 away, and relaxed from there the state runs away
    real = real_of(["s", "w"], {"p": 0.0, "al": 0.5}, {"s": "p*s - s**3 + 0.2*(w - al*p)", "w": "(w - al*p)*(-1 + 14*(w - al*p)**2)"}, "p")
    nf = an.normal_form(real, np.zeros(2), "p", 0.0)
    sweep = cm.Sweep(real, nf, RATE)
    q0, how = cm.start_of_sweep(real, sweep, np.zeros(2), cm.critical_space(real, nf, 1))
    assert how.startswith("the unwritten branch")
    assert np.allclose(q0, [0.0, 0.5 * sweep.value(0.0)], atol=1e-6)


def test_a_start_that_does_not_relax_is_reported():
    # the relaxing mode is unstable at the start of the sweep and runs away from any displacement, so neither the branch nor the relaxation gives a start
    real = real_of(["s", "w"], {"p": 0.0, "al": 0.5}, {"s": "p*s - s**3", "w": "-(p + 0.3)*(w - al*p) + (w - al*p)**3"}, "p")
    law = law_at_zero(real)
    assert law["status"].startswith("the start of the sweep does not relax")


def test_the_reduced_law_follows_a_branch_that_sets_off_after_its_first_point():
    # on the reduced law's grid the branch does not move over the first step, which gave version 4's guard no scale, so that it refused every later step
    real = real_of(["x"], {"d": 0.0, "p": 6.0}, {"x": "p/(1 + (d - x)**2) - x"}, "d", -2.0, 12.0)
    sweep_q = []
    q, q_prev = np.array([0.5]), None
    R = np.ones((1, 1))
    for d in np.linspace(0.0, 1.0, 11):
        q_next, ok = cm.branch_step(real, q, q_prev, {"d": float(d)}, (R, R), max_step=np.inf)
        assert ok
        q, q_prev = q_next, q
        sweep_q.append(float(q[0]))
    assert np.max(np.abs([real.F(np.array([[x]]), d=float(d))[0, 0] for x, d in zip(sweep_q, np.linspace(0.0, 1.0, 11))])) < 1e-8
