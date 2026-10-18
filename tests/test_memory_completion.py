"""The law at an obstruction (memory.completion): the completion that each obstruction names, the law derived from the
sweep, and its certificate by the full stochastic dynamics.

Statistical thresholds only (the probabilities are Monte Carlo estimates on a few thousand samples); no solver
residual is pinned."""
from __future__ import annotations

from math import pi
from pathlib import Path

import numpy as np
import pytest

from fieldbridge.memory import analysis as an, completion, construct, spec

EX = Path(__file__).resolve().parents[1] / "examples" / "memory"


def _law(name: str, **kw):
    real = spec.load(EX / f"{name}.json")
    rep = completion.complete(real, np.random.default_rng(3), floors=(1.0,), **kw)
    return rep["write_points"]


def test_completion_kind_follows_the_obstruction():
    gen = {"kind": "supercritical pitchfork", "critical_dimension": 1}
    assert completion.completion_kind(gen, {"realizes_generator": True}) == "none"
    assert completion.completion_kind({"critical_dimension": 2}, {"degenerate": True}) == "critical space"
    assert completion.completion_kind({"critical_dimension": 1}, {"degenerate_cubic": True}) == "leading order"
    assert completion.completion_kind({"critical_dimension": 1}, {"asymmetric": True}) == "second parameter"
    assert completion.completion_kind({"critical_dimension": 1}, {"subcritical": True}) == "distant branch"
    assert completion.completion_kind({"hopf": True}, {"hopf": True}) == "phase"


def test_diverged_trajectories_are_not_assigned_to_a_state():
    real = spec.load(EX / "pitchfork.json")
    stable = [np.array([1.0]), np.array([-1.0])]
    finals = np.array([[0.9], [-1.1], [np.nan], [np.inf]])
    lab = completion.assign(real, finals, stable)
    assert lab.tolist() == [0, 1, -1, -1]
    assert completion._probs(lab, 2) == [0.5, 0.5]


def test_the_rung_derives_the_swept_pitchfork_law_where_nothing_obstructs():
    """Control: the law derived from the sweep alone (no constant put in) is the transferred one, pi^(1/4), with no
    offset, and the full dynamics certifies it."""
    pts = _law("pitchfork", n_samples=3000, n_traj=3000)
    assert len(pts) == 1 and pts[0]["completion"] == "none"
    L = pts[0]["law"]
    assert L["status"] == "derived" and L["certified"]
    for src in ("law", "certificate"):
        x = L["one_direction"][src][0]
        assert abs(x["response_constant"] - pi ** 0.25) < 0.12
        assert abs(x["offset"]) < 0.1


def test_a_subcritical_write_point_gets_a_certified_law():
    """Two equal tubes: a positive cubic used to stop the derivation; the decision is still made in the linear
    stage, and the law of the completed model (with the distant branch) is certified."""
    pts = _law("tubes", n_samples=3000, n_traj=3000)
    sub = [p for p in pts if p["completion"] == "distant branch"]
    assert sub, [p["kind"] for p in pts]
    L = sub[0]["law"]
    assert L["status"] == "derived" and L["certified"]
    assert abs(L["one_direction"]["certificate"][0]["response_constant"] - pi ** 0.25) < 0.15


def test_two_critical_directions_give_a_law_over_all_written_states():
    """A square-symmetric two-component model has two critical directions; the law over its four written states,
    for biases in eight directions of the critical plane, agrees with the full dynamics."""
    pts = _law("two_component_square", n_samples=2000, n_traj=2000)
    deg = [p for p in pts if p["completion"] == "critical space"]
    assert deg, [p["kind"] for p in pts]
    L = deg[0]["law"]
    assert L["critical_dimension"] == 2 and len(L["stable_states"]) == 4 and len(L["rows"]) == 8
    for row in L["rows"]:
        assert abs(sum(row["law"]) - 1) < 1e-9
        assert row["deviation"] < 0.06
    # a bias favours the states on its side of the plane: the most probable state changes with the direction
    assert len({int(np.argmax(r["law"])) for r in L["rows"]}) == 4
