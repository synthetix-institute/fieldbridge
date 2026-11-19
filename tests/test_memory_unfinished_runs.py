"""A stable focus that turns much faster than it decays, and an implicit run of a relaxation that ends without settling
while it approaches a stable zero (correction 34 of the discovery engine).

Explicit Euler multiplies an eigenvalue l = a + i b by |1 + l dt|. At the step 0.3 / |l| of safe_dt a fast-turning
stable focus grew at every step, and the explicit prefix of the relaxation carried its starts to values that were not
finite. The step is now bounded by |a| / |l|^2 over the stable eigenvalues, which binds only on spirals.

A start has settled when its speed falls below 1e-5. A stable focus that turns at a rate w much faster than it decays
meets that criterion only at a distance of order 1e-5 / w from the state, after hundreds of thousands of turns, and its
implicit runs end at the budget of correction 30 unsettled: the card stored no state. A run that ends unsettled now counts
as settled at the zero Newton reaches from its end, when that zero has no growing direction, is not a center, and the run
moved towards it. A run along a cycle around an unstable zero, around a center, or towards no zero is not settled."""
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


# a linear stable focus that turns at 3e4 and decays at 0.3
FOCUS = _spec("fast-turning stable focus", ["x", "y"], {"w": 3e4}, {"x": "-0.3*x - w*y", "y": "w*x - 0.3*y"}, "w",
              [2e4, 4e4])
# a center that turns at 3e4: no state, and BDF damps it numerically
CENTER = _spec("fast-turning center", ["x", "y"], {"w": 3e4}, {"x": "-w*y", "y": "w*x"}, "w", [2e4, 4e4])
# a nonlinear center, x'' = -K x^3, whose origin has two flat rates: no state, and BDF damps it numerically
NONLINEAR_CENTER = _spec("fast nonlinear center", ["x", "y"], {"K": 1e9}, {"x": "y", "y": "-K*x**3"}, "K", [5e8, 2e9])
# a stable cycle of radius 1 around an unstable origin, turning at 3e4
CYCLE = _spec("fast-turning cycle", ["x", "y"], {"w": 3e4},
              {"x": "x*(1 - x**2 - y**2) - w*y", "y": "y*(1 - x**2 - y**2) + w*x"}, "w", [2e4, 4e4])


def _states(s, n=8):
    return an.stored_states(spec.load(s), np.random.default_rng(0), n)


def test_a_fast_turning_stable_focus_is_one_state():
    states, counts, spectra, unconverged, continuum = _states(FOCUS)
    assert len(states) == 1 and unconverged == 0
    assert np.linalg.norm(states[0]) < 1e-6


def test_a_fast_turning_center_stores_nothing():
    states, counts, spectra, unconverged, continuum = _states(CENTER)
    assert states == []


def test_a_fast_nonlinear_center_stores_nothing():
    real = spec.load(NONLINEAR_CENTER)
    assert 300.0 / an.safe_dt(real) > an.RELAX_BUDGET  # its relaxation is implicit
    states, counts, spectra, unconverged, continuum = _states(NONLINEAR_CENTER)
    assert states == []


def test_a_fast_cycle_around_an_unstable_zero_stores_nothing():
    states, counts, spectra, unconverged, continuum = _states(CYCLE)
    assert states == []


def test_the_step_keeps_a_fast_turning_focus_stable():
    # first part of the correction: Euler keeps l = a + i b decaying only for dt < 2 |a| / |l|^2
    real = spec.load(FOCUS)
    dt = an.safe_dt(real)
    w = 4e4  # the fastest rotation over the control range
    assert dt <= 0.3 / (0.09 + w ** 2) * (1 + 1e-6)
    assert abs(1 + complex(-0.3, w) * dt) < 1.0


def test_the_step_of_real_rates_is_unchanged():
    stiff = _spec("stiff relaxation", ["x", "y"], {"K": 1e6}, {"x": "-x + y", "y": "-K*(y - 0.5*x)"}, "K", [5e5, 2e6])
    real = spec.load(stiff)
    rates = [abs(np.linalg.eigvals(an.jacobian(real, np.zeros(2), {"K": k}))).max() for k in (5e5, 2e6)]
    assert abs(an.safe_dt(real) - 0.3 / max(rates)) < 1e-6 * an.safe_dt(real)
    plain = _spec("pitchfork", ["x"], {"r": 0.3}, {"x": "r*x - x**3"}, "r", [-0.5, 0.5])
    assert an.safe_dt(spec.load(plain)) == 0.01
