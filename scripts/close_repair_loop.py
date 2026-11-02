#!/usr/bin/env python
"""The first repair campaign on a records file: every law departure gets a task and one withheld test.

Retains all outcomes in a fresh library, runs the order closure on the laws of the selected write classes, and adds
three controls: the exact identity search on the same tasks, the same closure on the certified laws, and the number of
development evaluations with and without a hint from earlier verified campaigns of the class.

    python scripts/close_repair_loop.py docs/records/quantum_semiclassical.jsonl --out-dir OUT

--second-pass DIR runs again only the laws of a finished campaign that stopped before their withheld test, with the same
thresholds, and writes second_pass.json beside the first result, which it leaves unchanged.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import signal
import time

import numpy as np

from fieldbridge.memory import learning as ln
from fieldbridge.memory import repair_families as rf
from fieldbridge.memory import repair_loop as rl

CAMPAIGN = ("nonlinear stage", "equilibrium write")
TRIM = ("candidate", "policy", "record", "write_point", "failed_law")


def _budget(*_):
    raise TimeoutError("the time budget of this control ended")


def exact_baseline(task, library, seconds):
    """The exact identity on the same task, through the unchanged search and the direct-symbolic proposer."""
    out = {}
    for name, run in (("search", lambda: ln.learn(task, ln.MoveLibrary(library / "search.sqlite"),
                                                  ln.Policy(max_candidates=20))),
                      ("symbolic", lambda: rf.construct(task, None, ln.Policy(max_candidates=20),
                                                        evaluation_library=ln.MoveLibrary(library / "symbolic.sqlite"),
                                                        strategy="symbolic"))):
        started = time.monotonic()
        signal.signal(signal.SIGALRM, _budget)
        signal.alarm(int(seconds))
        try:
            report = run()
            out[name] = {"status": report["status"], "candidate_evaluations": report.get("candidate_evaluations"),
                         "reason": report.get("reason"), "description": rl.assess_exact_report(task, report),
                         "outcomes": dict(Counter(t.get("outcome", "") for t in report.get("trials", [])))}
        except TimeoutError as error:
            out[name] = {"status": "timed_out", "reason": str(error)}
        finally:
            signal.alarm(0)
        out[name]["seconds"] = round(time.monotonic() - started, 2)
    return out


def direct_evaluations(report):
    """Development checks when the ladder stops at the first order that passes."""
    rows = report.get("development_evaluations", [])
    used = 0
    for row in rows:
        used += 1
        if row["passed"]:
            break
    return used


def hinted_evaluations(report, hint):
    """Development checks when an earlier verified campaign of the class proposes its order and radius first; on a
    failure the ladder runs as without it. A pass of the hint can select another pair than the ladder does."""
    rows = report.get("development_evaluations", [])
    known = {(r["radius"], r["order"]): r["passed"] for r in rows}
    if hint is None:
        return direct_evaluations(report), True
    if known.get(hint):
        return 1, hint == (report.get("selected_radius"), report.get("selected_order"))
    return 1 + direct_evaluations(report), True


def _model_of(root, report):
    return {**json.loads((root / report["record"]["source"]).read_text(encoding="utf-8")), "noise": 0}


def second_pass(args):
    policy = rl.OrderPolicy(max_order=args.max_order)
    out = Path(args.second_pass)
    root = Path(args.root).resolve()
    records_path = Path(args.records)
    first = json.loads((out / "repair_loop.json").read_text(encoding="utf-8"))
    result = {"schema": "fieldbridge-repair-loop-second-pass/1", "policy": rl.asdict(policy),
              "implementation_sha256": rl.hashlib.sha256(Path(rl.__file__).read_bytes()).hexdigest(),
              "campaign": [], "certified_laws": [], "exact_baseline": [], "audit_of_first_pass": []}
    library = ln.MoveLibrary(out / "campaign.sqlite")
    # the verifier of the first pass scaled its limit to the first omitted coefficient: the selected candidates again
    for row in first["campaign"]:
        report = library.campaign(row["departure"])
        if report.get("candidate"):
            check = rl.verify_modulo_order({"model": _model_of(root, report)}, rl.candidate_from_record(report["candidate"]),
                                           policy)
            result["audit_of_first_pass"].append({"departure": row["departure"], "source": row["source"],
                                                  "target": row["target"], "status": row["status"], **check})
            print(json.dumps({"audit": row["source"], "target": row["target"], "first_status": row["status"],
                              "passed": check["passed"]}), flush=True)
    again = {r["departure"] for r in first["campaign"] if r["status"] == "unsupported_or_resonant"}
    if again:
        rows = rl.close_departures(library, root=root, policy=policy, select=set(args.classes),
                                   revisit=("unsupported_or_resonant",))
        for row in rows:
            if row["departure"] in again:
                report = library.campaign(row["departure"])
                result["campaign"].append({"departure": row["departure"],
                                           **{a: b for a, b in report.items() if a not in TRIM[:3]}})
                print(json.dumps({"campaign": report["failed_law"]["target"], "source": report["record"]["source"],
                                  "status": report["status"]}), flush=True)
    records = {r["source"]: r for r in map(json.loads, filter(str.strip, records_path.read_text(encoding="utf-8").splitlines()))}
    for row in first["campaign"]:
        statuses = {v.get("status") for v in (row.get("exact_baseline") or {}).values()}
        # the coordinate name of the generator collided with a variable of this model
        if "unsupported_or_invalid" in statuses and "non-reserved" in json.dumps(row["exact_baseline"]):
            task = rl.learning_task(records[row["source"]], row["target"], root=root, policy=policy)
            baseline = exact_baseline(task, _control_dir(out / "second", row["departure"]), args.baseline_seconds)
            result["exact_baseline"].append({"source": row["source"], "target": row["target"], **baseline})
            print(json.dumps({"exact": row["target"], "source": row["source"],
                              **{k: v["status"] for k, v in baseline.items()}}), flush=True)
    library2 = ln.MoveLibrary(out / "certified.sqlite")
    for row in first.get("certified_laws", []):
        if row["status"] == "unsupported_or_resonant" and "resonance" in str(row.get("reason")):
            task = rl.learning_task(records[row["record"]], row["target"], root=root, policy=policy)
            report = rl.close_modulo_order(task, library2, policy)
            result["certified_laws"].append({"record": row["record"], "target": row["target"],
                                             **{k: v for k, v in report.items() if k not in TRIM}})
            print(json.dumps({"certified": row["target"], "record": row["record"], "status": report["status"]}),
                  flush=True)
    text = json.dumps(rl._finite(result), indent=1, allow_nan=False) + "\n"
    (out / "second_pass.json").write_text(text, encoding="utf-8")


def run(args):
    policy = rl.OrderPolicy(max_order=args.max_order)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    root = Path(args.root).resolve()
    records_path = Path(args.records)
    records = [json.loads(line) for line in records_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    library = ln.MoveLibrary(out / "campaign.sqlite")
    retained = Counter()
    for record in records:
        retained.update(r["kind"] for r in ln.retain_record(library, record, source=str(records_path)))
    rows = rl.close_departures(library, root=root, policy=policy, select=set(args.classes),
                               law_seconds=args.law_seconds)
    reports = {}
    for row in rows:
        if row["status"] not in {"no_failed_law", "outside_this_campaign"}:
            reports[row["departure"]] = library.campaign(row["departure"])
    print(json.dumps({"retained": retained, "loop": rl.summary(rows)}), flush=True)

    # the exact identity on the same tasks
    tasks = {}
    if args.baseline_seconds > 0:
        for key, departure in library.departures():
            if key not in reports:
                continue
            record = rl._record_of(departure, root)
            try:
                task = rl.learning_task(record, departure["target"], root=root, policy=policy)
            except (ValueError, KeyError, ArithmeticError, np.linalg.LinAlgError, TypeError):
                continue
            tasks[key] = task
            reports[key]["exact_baseline"] = exact_baseline(task, _control_dir(out, key), args.baseline_seconds)
            print(json.dumps({"exact": departure["target"], "name": record["name"][:40],
                              **{k: v["status"] for k, v in reports[key]["exact_baseline"].items()}}), flush=True)

    # the certified laws, where nothing failed: does the closure separate the two groups?
    certified = []
    if args.certified:
        library2 = ln.MoveLibrary(out / "certified.sqlite")
        for record in records:
            for target, _, _, law in rl._locate_laws(record):
                if law.get("write_class") != "linear-stage write" or law.get("critical_dimension") != 1:
                    continue
                signal.signal(signal.SIGALRM, _budget)
                signal.alarm(args.law_seconds)
                try:
                    task = rl.learning_task(record, target, root=root, policy=policy)
                    report = rl.close_modulo_order(task, library2, policy)
                except TimeoutError as error:
                    report = {"status": "timed_out", "reason": str(error)}
                except (ValueError, KeyError, ArithmeticError, np.linalg.LinAlgError, TypeError) as error:
                    report = {"status": "no_task_generated", "reason": f"{type(error).__name__}: {error}"}
                finally:
                    signal.alarm(0)
                certified.append({"record": record["source"], "target": target,
                                  **{k: v for k, v in report.items() if k not in TRIM}})
                print(json.dumps({"certified": target, "name": record["name"][:40], "status": report["status"],
                                  "order": report.get("selected_order"), "radius": report.get("selected_radius")}),
                      flush=True)

    # memory against the direct ladder
    hints, memory = {}, []
    for key in sorted(reports):
        report = reports[key]
        cls = (report.get("failed_law") or {}).get("write_class")
        if "development_evaluations" in report:
            used, same = hinted_evaluations(report, hints.get(cls))
            memory.append({"departure": key, "class": cls, "direct": direct_evaluations(report), "hinted": used,
                           "same_selection": same, "hint": hints.get(cls)})
        if report.get("status") == "verified_smaller_description":
            hints[cls] = (report["selected_radius"], report["selected_order"])

    result = {"schema": "fieldbridge-repair-loop/1", "records": str(records_path), "policy": rl.asdict(policy),
              "implementation_sha256": rl.hashlib.sha256(Path(rl.__file__).read_bytes()).hexdigest(),
              "learning_sha256": rl.hashlib.sha256(Path(ln.__file__).read_bytes()).hexdigest(),
              "repair_families_sha256": rl.hashlib.sha256(Path(rf.__file__).read_bytes()).hexdigest(),
              "retained_outcomes": dict(retained), "total_retained": sum(retained.values()), "loop": rl.summary(rows),
              "campaign": [{"departure": k, **{a: b for a, b in v.items() if a not in TRIM},
                            "source": (v.get("record") or {}).get("source"),
                            "target": (v.get("failed_law") or {}).get("target"),
                            "write_class": (v.get("failed_law") or {}).get("write_class"),
                            "write_point": v.get("write_point")} for k, v in sorted(reports.items())],
              "certified_laws": certified, "memory_against_direct": memory, "novelty_established": False}
    text = json.dumps(rl._finite(result), indent=1, allow_nan=False) + "\n"
    (out / "repair_loop.json").write_text(text, encoding="utf-8")
    print(json.dumps({"status": Counter(r.get("status") for r in reports.values()),
                      "certified": Counter(r.get("status") for r in certified),
                      "output": str(out / "repair_loop.json")}), flush=True)


def _control_dir(out, key):
    path = out / "exact" / key[:16]
    path.mkdir(parents=True, exist_ok=True)
    return path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("records")
    parser.add_argument("--out-dir")
    parser.add_argument("--root", default=".")
    parser.add_argument("--classes", nargs="+", default=list(CAMPAIGN))
    parser.add_argument("--max-order", type=int, default=7)
    parser.add_argument("--baseline-seconds", type=int, default=120)
    parser.add_argument("--certified", action="store_true")
    parser.add_argument("--law-seconds", type=int, default=900, help="Budget of one law of the campaign and of the certified control")
    parser.add_argument("--second-pass", metavar="DIR")
    arguments = parser.parse_args()
    (second_pass if arguments.second_pass else run)(arguments)
