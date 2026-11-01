"""The linear stage of a completed write about the state the sweep really follows (completion.linear_stage_law).

The numbers come from the equations: exact Gaussian moments with Jacobians written by hand, computed outside this package, and closed
forms in the sweep rate. The law's moments carry no sampling error; the gate is the Euler step of their integration."""
from __future__ import annotations

from math import pi, sqrt

import numpy as np

from fieldbridge.memory import analysis as an, completion as cm, spec

RATE = 0.02
PI4 = pi ** 0.25
TOOL = {"T": 0.0, "B": 1.0, "G0": 0.5, "Ea": 3.75}
RELAX = "G0*exp(Ea*T)*(T - f)"


def real_of(variables, parameters, drift, control, lo=-0.8, hi=0.8, dt=0.02):
    return spec.load({"schema": "fieldbridge-memory/1", "kind": "equations", "name": "lagging state", "field": "test", "question": "what the write follows", "assumptions": [], "provenance": {"source": "test"},
                      "carrier": {"kind": "euclid", "name": "state", "variables": list(variables), "scale": 1.0}, "parameters": parameters,
                      "drift": drift, "noise": 0.001, "dt": dt, "control": {"name": control, "range": [lo, hi]}})


def stage(real, v=0.0, k=None):
    """The linear stage at the write point at v, for a bias of one floor toward each side of the first critical direction."""
    rng = np.random.default_rng(0)
    q = np.zeros(real.carrier.dim)
    events = an.locate_writes(real, rng, n_starts=12)
    ev = min(events, key=lambda e: abs(float(e["v"]) - v))
    nf = an.normal_form(real, ev["q"], real.control, ev["v"])
    k_static = max(1, int(nf.get("critical_dimension", 1)))
    k = k or max(k_static, an.dynamical_critical_dimension(nf, RATE, cm.SLOW_FACTOR))
    R, L = cm.critical_space(real, nf, k)
    sweep = cm.Sweep(real, nf, RATE)
    qc = np.asarray(nf["state"], float)
    q0 = cm._relaxed(real, qc, sweep.over(0.0))
    stable = cm.stable_states(real, rng, sweep.over(sweep.T), seeds=cm.swept_seeds(real, sweep, q0, R, sweep.T))
    D = 0.05 * RATE
    floor = sqrt(D) * RATE ** 0.25 / PI4
    axis = next(j for j in range(R.shape[0]) if real.variables[j] in ("s", "x"))      # the direction that crosses
    e = R.T @ np.eye(R.shape[0])[axis]
    e = e / np.linalg.norm(e)
    targets = [qc + 1e3 * sg * (R @ e) for sg in (1, -1)]
    crossing = cm.critical_space(real, nf, k_static) if k > k_static else None
    rows = cm.linear_stage_law(real, sweep, q0, R, L, targets, [floor, floor], D, stable, 200, rng, crossing=crossing)
    z = [float(e @ np.asarray(r["mean_critical"])) / sqrt(float(e @ np.asarray(r["cov_critical"]) @ e)) for r in rows]
    return {"c": PI4 * 0.5 * (z[0] - z[1]), "z0": 0.5 * (z[0] + z[1]), "t_map": rows[0]["t_map"], "t_cross": rows[0]["t_cross"],
            "lag": rows[0]["lag_at_crossing"], "rate": rows[0]["rate_at_crossing"], "k": k}


def test_a_lag_that_grows_along_the_sweep_slows_the_critical_rate_and_raises_the_constant():
    r = stage(real_of(["s", "f"], TOOL, {"s": "-B*f*s - s**3", "f": RELAX}, "T"))
    assert abs(r["c"] - 1.3826) < 0.01                       # first order in r: pi^(1/4) (1 - r E / (B G0))^(-1/4) = 1.387
    assert 0.8 < r["rate"] / RATE < 0.9                      # 1 - r E / (B G0) = 0.85
    assert r["lag"] > 0.02


def test_a_lag_that_shrinks_along_the_sweep_lowers_the_constant():
    r = stage(real_of(["s", "f"], TOOL, {"s": "B*f*s - s**3", "f": RELAX}, "T"))
    assert abs(r["c"] - 1.2737) < 0.01                       # first order in r: pi^(1/4) (1 + r E / (B G0))^(-1/4) = 1.286
    assert 1.1 < r["rate"] / RATE < 1.25


def test_a_state_that_cannot_lag_keeps_the_constant_and_the_declared_rate():
    r = stage(real_of(["s", "y"], {"p": 0.0, "g": 0.3}, {"s": "p*s - s**3 - s*y", "y": "-g*y + s**2"}, "p", -1.0, 1.0))
    assert abs(r["c"] - 1.3303) < 0.01
    assert r["lag"] < 1e-9
    assert abs(r["rate"] / RATE - 1) < 0.01


def test_a_constant_lag_delays_the_crossing_and_keeps_the_rate():
    g = 0.5
    r = stage(real_of(["s", "f"], {"T": 0.0, "B": 1.0, "g": g}, {"s": "-B*f*s - s**3", "f": "g*(T - f)"}, "T"))
    assert abs(r["c"] - 1.3303) < 0.01
    assert abs(r["t_map"] - r["t_cross"] - (1 / g + sqrt(2 / RATE))) < 0.2
    assert abs(r["lag"] - RATE / g) < 0.05 * RATE / g
    assert abs(r["rate"] / RATE - 1) < 0.01


def test_the_motion_of_the_branch_along_the_critical_direction_stays_a_forcing_on_the_mean():
    r = stage(real_of(["x"], {"p": 0.0, "alpha": 0.3}, {"x": "p*(x - alpha*p) - (x - alpha*p)**3"}, "p", -1.0, 1.0))
    assert abs(r["c"] - 1.3303) < 0.01
    assert abs(abs(r["z0"]) - 0.671) < 0.02
    assert abs(r["rate"] / RATE - 1) < 0.01


def test_a_slow_mode_inside_the_critical_space_follows_the_flow():
    r = stage(real_of(["s", "f"], {"T": 0.0, "B": 1.0, "G0": 0.1, "Ea": 0.6}, {"s": "-B*f*s - s**3", "f": RELAX}, "T", -0.9, 0.9))
    assert r["k"] == 2
    assert abs(r["c"] - 1.3640) < 0.01
    assert r["lag"] > 0.1                                    # about r / G0 = 0.2 at the crossing of the static branch
    assert r["t_map"] - r["t_cross"] > 15                    # the crossing comes about 1 / G0 later than the static one
