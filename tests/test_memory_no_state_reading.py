"""The reading of a card that stores no state (correction 29 of the discovery engine).

loss_law read "no stable state: a limit cycle, whose phase diffuses" whenever no state was stored and some starts did
not settle, whatever the trajectories did: a runaway, a blow-up in finite time or a chaotic motion read as a limit
cycle. The card now reads a limit cycle only where its search for a cycle finds one, neutral cycles where the Floquet
multiplier marks a family of them, and otherwise that no state is stored and the trajectories do not settle."""
from __future__ import annotations

import numpy as np

from fieldbridge.memory import discovery
from fieldbridge.memory import spec


def _spec(name, variables, parameters, drift, control, control_range, scale=1.0):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": name, "question": "what does it store?",
            "assumptions": ["test model"], "provenance": {"source": "test"},
            "carrier": {"kind": "euclid", "name": "state", "variables": variables, "scale": scale},
            "parameters": parameters, "drift": drift, "control": {"name": control, "range": control_range},
            "noise": 0.01, "dt": 0.01}


RUNAWAY = _spec("unstable node", ["x", "y"], {"a": 1.0}, {"x": "a*x", "y": "a*y + 0.5*x"}, "a", [0.5, 1.5])
BLOWUP = _spec("blow-up in finite time", ["x"], {"c": 1.0}, {"x": "c + x**2"}, "c", [0.5, 1.5])
VAN_DER_POL = _spec("van der Pol oscillator", ["x", "y"], {"mu": 1.0}, {"x": "y", "y": "mu*(1 - x**2)*y - x"}, "mu",
                    [0.5, 1.5], scale=2.0)
LOTKA_VOLTERRA = _spec("Lotka-Volterra", ["x", "y"], {"a": 1.0}, {"x": "x*(a - y)", "y": "y*(x - 1)"}, "a",
                       [0.5, 1.5], scale=2.0)
LOTKA_VOLTERRA["carrier"]["kind"] = "orthant"


def _card(s):
    return discovery.evaluate(spec.load(s), np.random.default_rng(0), quick=True)


def test_a_runaway_is_not_read_as_a_limit_cycle():
    card = _card(RUNAWAY)
    assert card["states"]["count"] == 0
    assert "limit cycle" not in card["loss_law"]
    assert discovery.mechanism_class(card) == "no stable state"


def test_a_blow_up_in_one_variable_is_not_read_as_a_limit_cycle():
    card = _card(BLOWUP)
    assert card["states"]["count"] == 0
    assert discovery.mechanism_class(card) == "no stable state"


def test_a_limit_cycle_that_the_search_finds_is_still_read():
    card = _card(VAN_DER_POL)
    assert card["loss_law"] == "no stable state: a limit cycle, whose phase diffuses (Law 2 along the cycle)"
    assert discovery.mechanism_class(card) == "limit cycle: phase memory"
    assert card["states"]["transverse_floquet_multiplier"] < 0.5


def test_neutral_cycles_are_still_read():
    card = _card(LOTKA_VOLTERRA)
    assert discovery.mechanism_class(card) == "neutral cycles: no isolated phase"
