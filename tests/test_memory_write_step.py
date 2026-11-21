"""The Euler-Maruyama step of the swept-write check (construct.write_step).

safe_dt bounds the Euler step at a weakly damped focus over random states of the carrier (correction 34 of the
discovery engine). The write check needs the bound only on the states its trajectories pass through: the critical
state along the sweep and the paths to the written states. A pitchfork in x beside a pair (u, v) that turns at the
rate 40 y^2 has a focus only away from y = 0, where the write runs."""
from __future__ import annotations

import numpy as np

from fieldbridge.memory import analysis as an
from fieldbridge.memory import construct as cs
from fieldbridge.memory import spec

FAR_FOCUS = {"schema": "fieldbridge-memory/1", "kind": "equations", "name": "pitchfork beside a far focus",
             "question": "test", "assumptions": ["test model"], "provenance": {"source": "test"},
             "carrier": {"kind": "euclid", "name": "state", "variables": ["x", "y", "u", "v"], "scale": 2.0},
             "parameters": {"eps": 0.0, "g": 0.5},
             "drift": {"x": "eps*x - x**3", "y": "-y", "u": "-g*u + 40*y**2*v", "v": "-40*y**2*u - g*v"},
             "control": {"name": "eps", "range": [-0.5, 0.5]}, "noise": 0.01, "dt": 0.01}


def _sweep(real):
    nf = an.normal_form(real, np.zeros(4), "eps", 0.0)
    rate = 0.02
    t_half = cs.sweep_half_time(real, nf, rate)
    v0, _ = cs.sweep_ends(nf, rate, t_half)
    dv = np.sign(-nf["kappa_slope"]) * rate / abs(nf["kappa_slope"])
    return nf, (lambda t: v0 + dv * t), t_half


def test_a_focus_away_from_the_write_leaves_its_step():
    real = spec.load(FAR_FOCUS)
    nf, path, t_half = _sweep(real)
    w = cs.write_step(real, nf, path, t_half)
    assert w == an.rate_dt(real)
    assert an.safe_dt(real) < 0.05 * w


def test_a_focus_on_the_write_bounds_its_step():
    """With the pair turning at 40 everywhere, the bound |a| / |l|^2 = 0.5 / (0.25 + 1600) holds on the write too."""
    real = spec.load(dict(FAR_FOCUS, drift=dict(FAR_FOCUS["drift"], u="-g*u + 40*v", v="-40*u - g*v")))
    nf, path, t_half = _sweep(real)
    w = cs.write_step(real, nf, path, t_half)
    assert abs(w - 0.5 / (0.25 + 1600.0)) < 1e-9


def test_refining_divides_the_write_step():
    real = spec.load(dict(FAR_FOCUS, drift=dict(FAR_FOCUS["drift"], u="-g*u + 40*v", v="-40*u - g*v")))
    nf, path, t_half = _sweep(real)
    w = cs.write_step(real, nf, path, t_half)
    an.refine_steps(real, 2.0)
    assert abs(cs.write_step(real, nf, path, t_half) - w / 2) < 1e-12


def test_a_focus_that_needs_too_short_a_step_is_not_run():
    """Turning at 400 with damping 0.5 the bound is 0.5 / 160000 = 3e-6, against 0.3 / 400 = 7.5e-4 for the rate: the
    check returns the reason and no measurement."""
    real = spec.load(dict(FAR_FOCUS, drift=dict(FAR_FOCUS["drift"], u="-g*u + 400*v", v="-400*u - g*v")))
    nf = an.normal_form(real, np.zeros(4), "eps", 0.0)
    assert cs.obstruction(nf)["realizes_generator"]
    law = cs.swept_write_check(real, nf, np.random.default_rng(0), n_traj=20)
    assert "measured" not in law and "focus" in law["unresolved"]
    assert law["rate_dt"] / law["dt"] > cs.MAX_WRITE_REFINEMENT
