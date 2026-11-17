"""Stored states that no barrier separates: a degenerate state is one state, and states on a set of equilibria are a
continuum.

The relaxation behind `stored_states` stops where the speed falls below its tolerance. Near a state whose restoring
force vanishes at linear order (x' = -x^3), the approach is a power of time, so starts from the two sides stop at
x = +-(1e-4)^(1/3) = +-0.046, farther apart than the merge tolerance, and were counted as two states. On a line of
equilibria (the epidemic of Kermack and McKendrick ends anywhere on i = 0) the starts merge into fewer clusters than
the count of a continuum requires. In both cases the card took two of these points as stored states, the write-field
continuation lost the first at zero field (it moves freely along its flat direction), and the card divided by that
zero field."""
from __future__ import annotations

import numpy as np
import pytest

from fieldbridge.memory import analysis as an
from fieldbridge.memory import discovery, spec


def _spec(name, variables, parameters, drift, control, control_range, carrier="euclid"):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": name, "question": "what does it store?",
            "assumptions": ["test model"], "provenance": {"source": "test"},
            "carrier": {"kind": carrier, "name": "state", "variables": variables, "scale": 1.0},
            "parameters": parameters, "drift": drift, "control": {"name": control, "range": control_range},
            "noise": 0.01, "dt": 0.01}


CUBIC = _spec("cubic relaxation with a slaved variable", ["x", "y"], {"g": 1.0}, {"x": "-g*x**3", "y": "-y + x**2"},
              "g", [0.5, 1.5])
SIR = _spec("Kermack-McKendrick epidemic", ["s", "i"], {"beta": 2.0, "gamma": 1.0},
            {"s": "-beta*s*i", "i": "beta*s*i - gamma*i"}, "beta", [1.0, 3.0], carrier="orthant")
PITCHFORK = _spec("pitchfork past onset", ["x", "y"], {"eps": 0.04}, {"x": "eps*x - x**3", "y": "-y"}, "eps",
                  [0.02, 0.08])


@pytest.mark.parametrize("seed", [0, 2, 3])
def test_the_cubic_relaxation_has_one_state_at_the_origin(seed):
    """Before the correction these seeds gave two states near x = +-0.04 (11 of the seeds 0-11 did)."""
    real = spec.load(CUBIC)
    states, counts, spectra, unconverged, continuum = an.stored_states(real, np.random.default_rng(seed), 48)
    assert len(states) == 1 and counts == [48] and unconverged == 0 and not continuum
    assert np.linalg.norm(states[0]) < 1e-3
    # degenerate: the restoring rate vanishes at linear order
    assert an.classify(spectra[0])["flat"] == 1


def test_the_card_of_the_cubic_relaxation_counts_one_state_and_does_not_divide_by_zero():
    card = discovery.evaluate(spec.load(CUBIC), np.random.default_rng(0), quick=True)
    assert card["states"]["count"] == 1 and "writes" not in card
    assert discovery.mechanism_class(card) == "single state: relaxation"


@pytest.mark.parametrize("seed", [0, 1, 5])
def test_the_end_states_of_an_epidemic_are_a_continuum(seed):
    """Every s < gamma / beta with i = 0 is a stable state; the starts end on that line, 3 to 5 clusters of 48. Its
    restoring rate gamma - beta s vanishes at the threshold, so the relaxation stops at i ~ 1e-5 near it."""
    real = spec.load(SIR)
    states, _, _, _, continuum = an.stored_states(real, np.random.default_rng(seed), 48)
    assert len(states) >= 2 and continuum
    assert all(abs(q[1]) < 1e-3 and q[0] < real.params["gamma"] / real.params["beta"] for q in states)
    card = discovery.evaluate(real, np.random.default_rng(seed), quick=True)
    assert card["states"]["continuum"] and "writes" not in card
    assert discovery.mechanism_class(card) == "zero mode: diffusive retention"


def test_two_states_of_a_pitchfork_past_onset_stay_two_relaxed_states():
    """States that are not degenerate are kept as relaxed: the correction leaves them, and their card, unchanged."""
    real = spec.load(PITCHFORK)
    states, counts, spectra, _, continuum = an.stored_states(real, np.random.default_rng(0), 48)
    assert len(states) == 2 and sum(counts) == 48 and not continuum
    assert sorted(round(float(q[0]), 2) for q in states) == [-0.2, 0.2]
    kept, kept_counts = an._one_state_per_zero(real, states, counts, {})
    assert kept_counts == counts and all(a is b for a, b in zip(kept, states))


def test_an_isolated_degenerate_zero_is_not_a_set_of_equilibria():
    real = spec.load(CUBIC)
    assert not an._on_a_set_of_equilibria(real, np.array([0.04, 0.0016]), {})
    sir = spec.load(SIR)
    assert an._on_a_set_of_equilibria(sir, np.array([0.3, 0.0]), {})
