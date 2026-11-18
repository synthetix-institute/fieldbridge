"""A start that moves slowly far from any zero of the drift has not settled (correction 24 of the discovery engine).

The relaxation behind stored_states counts a start as settled when its speed falls below 1e-4. On the slow manifold of
a relaxation oscillation, or under a slow drift that no zero ends, the speed is that small far from any zero, and the
starts were read as stored states: 22 of them on the slow manifold of the van der Pol oscillator at mu = 1e5, whose
noisy hold time then ran at the step of the fast rate. A start now counts as settled when one step of Newton from it is
shorter than the merge tolerance and accounts for the drift, or else when Newton takes it to a zero within that
tolerance, or to a zero without a growing direction, which the start approaches slowly."""
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


# a fast relaxation onto the line x = y, along which y drifts at the constant speed s: no zero anywhere
DRIFT = _spec("slow drift along a line", ["x", "y"], {"K": 10.0, "s": 1e-5}, {"x": "-K*(x - y)", "y": "s"}, "s",
              [5e-6, 2e-5])
# FitzHugh-Nagumo with the time scales 1e5 apart: with Iapp = 0.5 a relaxation oscillation, with Iapp = 0 one excitable state
FHN = {"v": "K*(v - v**3/3 - w + Iapp)", "w": "(v + 0.7 - 0.8*w)/K"}
OSCILLATING = _spec("relaxation oscillation of FitzHugh-Nagumo", ["v", "w"], {"K": 1e5, "Iapp": 0.5}, FHN, "Iapp",
                    [0.4, 0.6], scale=2.0)
EXCITABLE = _spec("excitable state of FitzHugh-Nagumo", ["v", "w"], {"K": 1e5, "Iapp": 0.0}, FHN, "Iapp", [-0.1, 0.1],
                  scale=2.0)
DEGENERATE = _spec("degenerate state", ["x"], {"a": 1.0}, {"x": "-a*x**3"}, "a", [0.5, 1.5])


def _states(s):
    return an.stored_states(spec.load(s), np.random.default_rng(0), 48)


def test_a_drift_that_no_zero_ends_stores_nothing():
    states, counts, spectra, unconverged, continuum = _states(DRIFT)
    assert states == [] and unconverged == 48


def test_the_slow_manifold_of_a_relaxation_oscillation_stores_nothing():
    states, counts, spectra, unconverged, continuum = _states(OSCILLATING)
    assert states == [] and unconverged == 48


def test_a_slow_approach_to_a_stable_zero_is_one_state():
    # the starts reach the slow manifold and creep towards the rest state at the rate 1 / K; each is taken to that
    # zero, which has no growing direction, so they are one state there
    states, counts, spectra, unconverged, continuum = _states(EXCITABLE)
    assert len(states) == 1 and unconverged == 0
    v, w = states[0]
    assert abs(v - v**3 / 3 - w) < 1e-6 and abs(v + 0.7 - 0.8 * w) < 1e-6


def test_a_degenerate_state_is_still_one_state():
    states, counts, spectra, unconverged, continuum = _states(DEGENERATE)
    assert len(states) == 1 and unconverged == 0
