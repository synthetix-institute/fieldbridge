"""Typed moves on a drift specification: the term handling, the boundedness test, the candidates, the obstruction
classes, the closure test of compositions and the command."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from fieldbridge.memory import moves as mv

ROOT = Path(__file__).resolve().parents[1]
PITCHFORK = json.loads((ROOT / "examples" / "memory" / "pitchfork.json").read_text(encoding="utf-8"))


def plain_pitchfork(drift="eps*x - x**3", control=("eps", [-0.5, 1.5]), **params):
    spec = dict(PITCHFORK)
    spec.pop("potential", None)
    spec["drift"] = {"x": drift}
    spec["parameters"] = dict(PITCHFORK["parameters"], **params)
    spec["control"] = {"name": control[0], "range": list(control[1])}
    spec["material"] = sorted(set(spec["parameters"]) - {control[0]})
    return spec


def test_terms_split_at_top_level_signs_only_and_join_back():
    assert mv.split_terms("eps*x - x**3 + h") == ["eps*x", "- x**3", "+ h"]
    assert mv.split_terms("-(K*(x**2 + y**2) - p + Delta)*x - 0.5*kappa*y") == ["-(K*(x**2 + y**2) - p + Delta)*x", "- 0.5*kappa*y"]
    assert mv.split_terms("2e-3*x + exp(-x)*y") == ["2e-3*x", "+ exp(-x)*y"]
    terms = mv.split_terms("a*x - b*x**3 + c")
    assert mv.join_terms(terms) == "a*x - b*x**3 + c" and mv.join_terms(terms[1:]) == "- b*x**3 + c"


def test_obstruction_classes_name_what_a_completion_keeps():
    assert mv.obstruction_class("one-sided write (fold) and no material parameter reaches a cusp")["kind"] == "second parameter"
    assert mv.obstruction_class("positive cubic term: the write is symmetric but the state jumps")["kind"] == "distant branch"
    assert mv.obstruction_class("degenerate critical mode: 2 directions become critical together")["kind"] == "critical space"
    assert mv.obstruction_class("vanishing cubic term")["kind"] == "leading order"
    assert mv.obstruction_class("no stable state (the realization oscillates)")["kind"] == "phase"
    assert mv.obstruction_class("a single stable state: there is no second one")["kind"] == "none"
    assert mv.obstruction_class("something else entirely") == {"kind": "unclassified", "keeps": "something else entirely", "has_law": False}
    assert mv.obstruction_class("time limit of 60 s")["kind"] == "not derived"


def test_the_boundedness_test_keeps_the_pitchfork_and_rejects_it_without_its_cubic():
    ok, why = mv.bounded(plain_pitchfork())
    assert ok and why == ""
    ok, why = mv.bounded(plain_pitchfork("eps*x"))
    assert not ok and "unbounded" in why and "eps = 1.5" in why  # the growth is at the positive end of the range
    ok, why = mv.bounded(plain_pitchfork("eps*x + x**3"))        # the subcritical pitchfork without its saturation
    assert not ok


def test_the_candidates_are_the_four_moves_and_the_unbounded_ones_are_rejected():
    spec = plain_pitchfork("eps*x - x**3 + h", h=0.1)
    rejected = []
    cands = mv.candidates(spec, rejected=rejected)
    moves = {(c["move"], c["what"]) for c in cands}
    assert ("add_protocol_order", "control = h over [0.025, 0.4]") in moves
    assert ("modify_protocol_order", "control eps range wider: [-20.5, 21.5]") in moves
    assert ("project_operator_apparatus", "drop '+ h' from dx/dt") in moves
    assert ("project_operator_apparatus", "drop 'eps*x' from dx/dt") in moves
    assert all(c["move"] != "add_operator_apparatus" for c in cands)   # no donor
    dropped_cubic = [r for r in rejected if r["what"] == "drop '- x**3' from dx/dt"]
    assert dropped_cubic and "unbounded" in dropped_cubic[0]["why"]
    donor = {"name": "quintic donor", "drift": {"y": "a*y - b*y**5 + c*sin(y)"}, "parameters": {"a": 1.0, "b": 1.0, "c": 0.5}}
    with_donor = [c for c in mv.candidates(spec, donor, check_bounded=False) if c["move"] == "add_operator_apparatus"]
    whats = {c["what"] for c in with_donor}
    assert "add g_add*(1*x) to dx/dt (from quintic donor)" in whats
    assert "add g_add*(1*x**5) to dx/dt (from quintic donor)" in whats
    assert "add g_add*(1*sin(x)) to dx/dt (from quintic donor)" in whats
    added = next(c for c in with_donor if "x**5" in c["what"])
    assert added["spec"]["parameters"]["g_add"] == 1.0 and "g_add" in added["spec"]["material"]


def test_a_composition_closes_on_known_outcomes_and_opens_on_a_new_target_or_class():
    base = (("threshold-write",), (("symmetric-write", "second parameter"),))
    step = (("symmetric-write", "threshold-write"), ())
    assert mv.closes(step, [base, step]) == (True, "equals a known outcome")
    assert mv.closes((("symmetric-write",), ()), [base]) == (False, "reaches symmetric-write that no single step reaches")
    assert mv.closes((("threshold-write",), (("symmetric-write", "critical space"),)), [base]) == (
        False, "new obstruction class symmetric-write: critical space")
    assert mv.closes((("threshold-write",), (("symmetric-write", "not derived"),)), [base])[0]
    assert mv.closes(((), (("symmetric-write", "second parameter"), ("threshold-write", "none"))),
                     [base, ((), (("threshold-write", "none"),))]) == (True, "a combination of known classes")


def test_the_judge_reads_the_pitchfork_and_the_bias_as_the_derivations_do():
    der = mv.judge(plain_pitchfork())
    assert der["symmetric-write"]["status"].startswith("reached") and der["threshold-write"]["status"].startswith("reached")
    reached, classes = mv.summary_of(der)
    assert reached == {"symmetric-write", "threshold-write"} and classes == {}
    der = mv.judge(plain_pitchfork("eps*x + x**3 - x**5", control=("eps", [-0.5, 0.5])))
    assert der["symmetric-write"]["status"] == "obstructed"
    assert mv.summary_of(der)[1]["symmetric-write"] == "second parameter"   # a fold along the control


def test_the_moves_command_finds_the_step_that_removes_a_bias(tmp_path):
    """The biased pitchfork reaches the symmetric write after unfolding; dropping the bias term makes it reach it
    directly, which the derivation reports with a different word, so the step list stays empty while the judged list
    holds every move. The command's report carries both."""
    spec = plain_pitchfork("eps*x - x**3 + h", h=0.1)
    path = tmp_path / "biased.json"
    path.write_text(json.dumps(spec), encoding="utf-8")
    out = tmp_path / "moves"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "moves", str(path), "--out-dir", str(out), "--no-isolate",
           "--max-candidates", "4"]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    assert res.returncode == 0, res.stderr
    summary = json.loads(res.stdout)
    report = json.loads((out / "moves.json").read_text(encoding="utf-8"))
    assert report["command"] == "moves" and report["candidates"] == summary["candidates"] >= 4
    assert len(report["judged"]) == 4 and all(j["move"] in mv.MOVES for j in report["judged"])
    assert any(r["what"] == "drop '- x**3' from dx/dt" for r in report["rejected"])
    md = (out / "moves.md").read_text(encoding="utf-8")
    assert md.startswith("# Typed moves on") and "rejected as unbounded" in md
