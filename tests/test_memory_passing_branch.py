"""A start near a zero with a growing direction, whose own rates have none, has not settled (correction 26 of the
discovery engine).

The slow branch of a relaxation oscillation can pass within the merge tolerance of its unstable steady state. Starts on
that branch move slowly, below the speed tolerance, and one step of Newton from them is short: under correction 24 they
counted as settled near a zero, and since their own rates are those of the attracting branch, the card read them as a
stable state with a flat direction. The zero they are near grows. Such a start has not settled; a start whose own rates
show the growth, as one that stopped near a saddle, is settled as before and is not kept as a state."""
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


# FitzHugh-Nagumo with the time scales 1e5 apart. With Iapp = 0.33 the rest point v = -0.968 lies just past the lower
# knee (v = -1), on the repelling middle branch: a relaxation oscillation whose left branch, near the knee, lies within
# the merge tolerance (0.14 at the scale 2) of that unstable rest point. With Iapp = 0.25 the rest point lies
# just before the knee, on the attracting branch: one excitable state, approached along the same branch.
FHN = {"v": "K*(v - v**3/3 - w + Iapp)", "w": "(v + 0.7 - 0.8*w)/K"}
KNEE = _spec("relaxation oscillation past the knee", ["v", "w"], {"K": 1e5, "Iapp": 0.33}, FHN, "Iapp", [0.3, 0.36],
             scale=2.0)
BEFORE_KNEE = _spec("excitable state before the knee", ["v", "w"], {"K": 1e5, "Iapp": 0.25}, FHN, "Iapp", [0.2, 0.28],
                    scale=2.0)
# the two-variable Oregonator in the time unit of its fast variable, at f = 1.5, q = 2e-3: a relaxation oscillation whose
# left branch, x ~ q, passes within 0.01 to 0.07 of the unstable steady state x = z = 0.0098 for z below 0.07
OREGONATOR = _spec("Oregonator with a slow branch beside its unstable steady state", ["x", "z"],
                   {"f": 1.5, "q": 2e-3, "eps": 1e-5},
                   {"x": "x*(1 - x) - f*z*(x - q)/(x + q)", "z": "eps*(x - z)"}, "f", [1.2, 1.8])
OREGONATOR["carrier"]["kind"] = "orthant"
SADDLE = _spec("double well with a relaxation", ["x", "y"], {"a": 1.0}, {"x": "a*x - x**3", "y": "-y"}, "a", [0.5, 1.5])


def _states(s):
    return an.stored_states(spec.load(s), np.random.default_rng(0), 48)


def test_the_oscillation_past_the_knee_stores_nothing():
    # under either code: the random starts land on the two branches far from the knee
    states, counts, spectra, unconverged, continuum = _states(KNEE)
    assert states == [] and unconverged == 48


def test_the_slow_branch_of_the_oregonator_stores_nothing():
    # under correction 24 five starts on the left branch formed a stable state at (0.0023, 0.0195)
    states, counts, spectra, unconverged, continuum = _states(OREGONATOR)
    assert states == [] and unconverged == 48


def test_the_starts_on_that_branch_are_near_the_unstable_zero():
    # the case of the correction: a start on the left branch, 0.08 from the rest point, whose own rates have no growing
    # direction, and whose Newton step is short
    real = spec.load(KNEE)
    v = -1.05
    q = np.array([v, v - v**3 / 3 + 0.33])
    J = an.jacobian(real, q)
    assert an.classify(-np.linalg.eigvals(J))["unstable"] == 0
    z, ok = an.newton(real, q)
    assert ok and real.carrier.distance(z, q) < an._merge_tol(real)
    assert an.classify(an.kappa_spectrum(real, z))["unstable"] > 0
    assert not an._settled(real, q, None)


def test_an_excitable_state_on_the_same_branch_is_one_state():
    states, counts, spectra, unconverged, continuum = _states(BEFORE_KNEE)
    assert len(states) == 1
    v, w = states[0]
    assert abs(v - v**3 / 3 - w + 0.25) < 1e-6 and abs(v + 0.7 - 0.8 * w) < 1e-6


def test_a_start_whose_own_rates_grow_near_a_saddle_is_settled_as_before():
    real = spec.load(SADDLE)
    q = np.array([1e-7, 1e-7])
    assert an.classify(-np.linalg.eigvals(an.jacobian(real, q)))["unstable"] == 1
    assert an._settled(real, q, None)


def test_stable_states_are_unchanged():
    states, counts, spectra, unconverged, continuum = _states(SADDLE)
    assert len(states) == 2 and unconverged == 0
    assert sorted(np.sign(float(s[0])) for s in states) == [-1.0, 1.0]
    assert all(abs(abs(float(s[0])) - 1.0) < 1e-3 for s in states)
