#!/usr/bin/env python3
"""Learn repair families on one family of equations, freeze them, and test construction on another family.

Training family (normal-form amplitude equations, authored): higher-order saturation, a hidden exchange imbalance,
a hidden squared amplitude and a moving cubic branch. Test family (kinetics, mechanics and population dynamics,
authored, with different dimensions, degrees and variable names): a three-stage linear cascade observed at its first
stage, a three-body product feeding an observed species, an oscillator whose spring anchor moves, a population with
quartic crowding, and a designed counterexample whose omitted variable generates a tower that does not close.

Four constructors are compared on the test family, each with its memory frozen after training:
  blind      the enumeration of learning.py with an empty library;
  recipes    level 1, the learning.py move library of verified supports;
  symbolic   direct symbolic closure, without learned memory;
  families   level 2, the repair-family library.
Afterwards a separate revision step records the counterexample and narrows the family's assumptions.
No corpus sweep and no external model call.
"""
import argparse
import copy
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fieldbridge.memory import learning as ln  # noqa: E402
from fieldbridge.memory import repair_families as rf  # noqa: E402


def task(name, variables, params, drift, coords, cand, dev, wh, bounds=None, field="test family"):
    model = {"schema": "fieldbridge-memory/1", "name": name, "kind": "equations", "question": name,
             "provenance": {"kind": "authored verification example", "source": field},
             "assumptions": ["deterministic", "constant parameters during each intervention"],
             "carrier": {"kind": "euclid", "variables": variables, "scale": 1},
             "parameters": params, "drift": drift, "noise": 0}
    if bounds:
        model["bounds"] = bounds
    return {"schema": ln.SCHEMA, "model": model, "coordinates": coords, "candidate_drift": cand,
            "development": dev, "withheld": wh}


def case(i, init, params, T=3.0):
    return {"id": i, "initial": init, "parameters": params, "duration": T}


def training_tasks():
    ex = ROOT / "examples" / "learning"
    quintic = json.loads((ex / "quintic_saturation.json").read_text())
    hidden = json.loads((ex / "hidden_exchange.json").read_text())
    square = copy.deepcopy(hidden)
    square["model"]["drift"]["x"] = "-a*x+b*y**2"
    moving = json.loads((ex / "moving_branch.json").read_text())
    return {"quintic_saturation": quintic, "hidden_exchange": hidden, "hidden_square": square,
            "moving_branch": moving}


def test_tasks():
    return {
        "kinetic_cascade": task(
            "three-stage linear cascade observed at its first stage", ["A", "B", "C"],
            {"k1": 1.0, "k2": 0.7, "k3": 0.6, "k4": 0.5, "k5": 0.4},
            {"A": "-k1*A+k2*B", "B": "-k3*B+k4*C", "C": "-k5*C"}, {"s": "A"}, {"s": "-k1*s"},
            [case("dev-1", [0.1, 0.8, 1.0], {}), case("dev-2", [0.3, -0.5, 0.7], {"k2": 0.9})],
            [case("wh-1", [0.0, 1.2, -0.6], {"k4": 0.8, "k5": 0.3}, 4), case("wh-2", [0.2, 0.4, 1.3], {"k3": 1.1}, 4)],
            field="chemical kinetics"),
        "three_body_product": task(
            "observed species fed by a product of two decaying species", ["c1", "c2", "c3"],
            {"a": 1.0, "b": 0.8, "c": 0.5, "d": 0.3},
            {"c1": "-a*c1+b*c2*c3", "c2": "-c*c2", "c3": "-d*c3"}, {"n": "c1"}, {"n": "-a*n"},
            [case("dev-1", [0.1, 1.0, 1.2], {}), case("dev-2", [0.4, -0.8, 0.9], {"b": 1.1})],
            [case("wh-1", [0.0, 1.3, 1.1], {"c": 0.7, "d": 0.2}, 4), case("wh-2", [-0.2, 0.9, -1.0], {"a": 1.4}, 4)],
            field="reaction kinetics"),
        "moving_anchor_oscillator": task(
            "oscillator whose spring anchor moves at constant speed", ["q", "p", "s"],
            {"k": 2.0, "g": 0.1, "w": 0.3},
            {"q": "p", "p": "-k*(q-s)-g*p", "s": "w"}, {"e": "q-s", "m": "p"},
            {"e": "m", "m": "-k*e-g*m"},
            [case("dev-1", [0.5, 0.0, 0.0], {}), case("dev-2", [0.0, 0.4, 1.0], {"w": -0.2})],
            [case("wh-1", [0.3, -0.2, 0.5], {"w": 0.5, "k": 1.5}, 4), case("wh-2", [-0.4, 0.1, 0.0], {"w": -0.4}, 4)],
            field="mechanics"),
        "quartic_crowding": task(
            "population limited by four-body crowding", ["pop"], {"r": 0.6, "s": 1.0},
            {"pop": "r*pop-s*pop**4"}, {"n": "pop"}, {"n": "r*n-s*n**2"},
            [case("dev-1", [0.3], {}), case("dev-2", [0.6], {"r": 0.4, "s": 0.8})],
            [case("wh-1", [0.2], {"r": 0.9, "s": 1.3}, 4), case("wh-2", [0.8], {"r": 0.5, "s": 1.7}, 4)],
            bounds={"s": [0.1, 3]}, field="population dynamics"),
        "nonlinear_hidden_feedback": task(
            "observed amplitude fed by a mode that is itself driven by its square", ["x", "y"],
            {"a": 1.0, "b": 0.5, "c": 0.8},
            {"x": "-a*x+b*y**2", "y": "-c*y+x**2"}, {"u": "x"}, {"u": "-a*u"},
            [case("dev-1", [0.2, 0.5], {}), case("dev-2", [-0.3, 0.4], {"b": 0.7})],
            [case("wh-1", [0.1, 0.6], {"c": 1.1}, 4), case("wh-2", [0.3, -0.5], {"a": 1.3}, 4)],
            field="designed counterexample"),
        "mixed_hidden": task(
            "mixed retained and hidden interaction", ["q", "z"], {"a": 1.0, "b": .5},
            {"q": "-a*q+b*q*z", "z": "-z"}, {"obs": "q"}, {"obs": "-a*obs"},
            [case("d1", [.2, .4], {}), case("d2", [.3, -.2], {"b": .8})],
            [case("w1", [.5, .6], {"a": .7}), case("w2", [.4, -.3], {"a": 1.3})],
            field="multiplicative feedback"),
        "budget_limited_cascade": task(
            "four-stage cascade requires three hidden coordinates", ["q1", "q2", "q3", "q4"],
            {"a": 1.0, "b": .5},
            {"q1": "-a*q1+b*q2", "q2": "-a*q2+b*q3", "q3": "-a*q3+b*q4", "q4": "-a*q4"},
            {"obs": "q1"}, {"obs": "-a*obs"},
            [case("d1", [.1, .2, .3, .4], {}), case("d2", [.3, .1, -.2, .5], {"b": .8})],
            [case("w1", [.2, .3, -.4, .6], {"a": .7}), case("w2", [.4, -.2, .5, .1], {"a": 1.3})],
            field="explicit coordinate-budget control"),
    }


def row(r):
    return {"status": r["status"], "candidate_evaluations": r["candidate_evaluations"],
            "selection_origin": r.get("selection_origin"), "groebner_reductions": r.get("groebner_reductions"),
            "added": (r.get("candidate") or {}).get("coordinates"), "drift": (r.get("candidate") or {}).get("drift"),
            "withheld_error": (r.get("withheld") or {}).get("max_error"),
            "omission_error": (r.get("withheld") or {}).get("max_omission_error")}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    if out.exists() and any(out.iterdir()):
        raise SystemExit("Choose a fresh directory")
    out.mkdir(parents=True, exist_ok=True)

    # ---- training: blind learning with storage of level-1 recipes, then abstraction into level-2 families
    recipes = ln.MoveLibrary(out / "train_moves.sqlite")
    families = rf.FamilyLibrary(out / "families.json")
    training = {}
    for name, t in training_tasks().items():
        rep = ln.learn(t, recipes)
        training[name] = {**row(rep), "family": families.abstract(t, rep, name),
                          "classification": rf.classify(t)}
        if rep["status"] != "verified_extension":
            raise SystemExit(f"training case {name} not verified: {rep['status']}")

    # A development failure revises hint selection, without reading any withheld outcomes.
    revision_task = training_tasks()["hidden_exchange"]
    limits = ln.Policy(max_added_coordinates=0)
    revision_before = rf._proposals(revision_task, families, limits)[0]
    development_failure = rf.construct(revision_task, families, limits,
                                       evaluation_library=ln.MoveLibrary(out / "revision.sqlite"))
    assert development_failure["status"] != "verified_extension"
    assert "withheld" not in development_failure
    families.revise("development-budget-control", revision_task, "omitted_state_variable",
                    "this task did not close within the declared coordinate budget", policy=limits)
    revision_after, blocked = rf._proposals(revision_task, families, limits)
    assert not blocked and "guidance_skipped" in revision_after
    revisions = [{"task": "development-budget-control", "before": revision_before,
                  "after": revision_after, "withheld_evaluated": False}]

    # ---- freeze both memories; each strategy has a separate, persistent evaluation registry
    frozen_digest = families.digest()
    frozen = rf.FamilyLibrary(out / "families.json", frozen=True)
    test_rows = {}
    for name, t in test_tasks().items():
        shutil.copy(out / "train_moves.sqlite", out / f"frozen_moves_{name}.sqlite")
        blind = ln.learn(t, ln.MoveLibrary(out / "blind_evaluations.sqlite"))
        level1 = ln.learn(t, ln.MoveLibrary(out / f"frozen_moves_{name}.sqlite"))
        direct = rf.construct(t, strategy="symbolic",
                              evaluation_library=ln.MoveLibrary(out / "symbolic_evaluations.sqlite"))
        level2 = rf.construct(t, frozen, evaluation_library=ln.MoveLibrary(out / "family_evaluations.sqlite"))
        test_rows[name] = {"classification": rf.classify(t), "blind": row(blind), "recipes": row(level1),
                           "symbolic": row(direct), "families": row(level2),
                           "family_analysis": level2.get("family_analysis")}
        for method, result in (("blind", blind), ("recipes", level1), ("symbolic", direct), ("families", level2)):
            test_rows[name][method]["seconds"] = result["seconds"]
    assert frozen.digest() == frozen_digest, "the frozen library changed during testing"

    def tally(key):
        rows = [r[key] for r in test_rows.values()]
        return {"verified": sum(r["status"] == "verified_extension" for r in rows),
                "candidate_evaluations": sum(r["candidate_evaluations"] for r in rows),
                "seconds": sum(r["seconds"] for r in rows),
                "groebner_reductions": sum(r.get("groebner_reductions") or 0 for r in rows),
                "tasks": len(rows)}

    summary = {"training": {k: {kk: v[kk] for kk in ("status", "candidate_evaluations", "family")}
                            for k, v in training.items()},
               "frozen_family_library_sha256": frozen_digest,
               "test": {k: {m: {kk: v[m][kk] for kk in ("status", "candidate_evaluations", "selection_origin")}
                            for m in ("blind", "recipes", "symbolic", "families")} for k, v in test_rows.items()},
               "totals": {m: tally(m) for m in ("blind", "recipes", "symbolic", "families")},
               "expected_unresolved": ["budget_limited_cascade"],
               "revisions": revisions}
    report = {"scope": "authored equation-level cases; tests transfer of repair families between equation families, "
                       "not novelty of any physical law", "training": training, "test": test_rows,
              "summary": summary}
    (out / "transfer_experiment.json").write_text(json.dumps(report, indent=2, allow_nan=False, default=str) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
