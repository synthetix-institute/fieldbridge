"""The record of a drift specification (memory.record, memory record): its steps, its assembly, the command."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from fieldbridge.memory import record as rc

ROOT = Path(__file__).resolve().parents[1]
PITCHFORK = json.loads((ROOT / "examples" / "memory" / "pitchfork.json").read_text(encoding="utf-8"))
FAST = {"law": False, "trajectories": 100, "samples": 400, "trajectories_completion": 300}


def test_the_write_points_step_reads_the_pitchfork():
    wp = rc.step(PITCHFORK, "write_points")
    assert wp["status"] == "derived" and wp["symmetry_order"] == 2 and wp["stable_states_at_declared_parameters"] == 2
    assert wp["control"] == {"name": "eps", "range": [-0.5, 1.5]}
    assert len(wp["points"]) == 1
    p = wp["points"][0]
    assert p["kind"].startswith("supercritical pitchfork") and p["critical_dimension"] == 1 and p["completion"] == "none"
    assert abs(p["control_value"]) < 1e-3 and abs(p["coefficients"]["a3"] + 1.0) < 0.05


def test_a_specification_without_a_control_is_recorded_as_such():
    s = dict(PITCHFORK)
    s.pop("control")
    s["material"] = []
    wp = rc.step(s, "write_points")
    assert wp["control"] is None and "no control parameter" in wp["note"] and wp["points"] == []
    assert rc.step(s, "completion")["status"] == "no control parameter"
    rec = rc.assemble(s, {"write_points": wp})
    assert rec["control"] is None and "no control parameter" in rec["note"]
    assert rec["targets"]["symmetric-write"] == {"status": "not run"}
    assert "reaches no write target" in rec["summary"]


def test_the_whole_record_of_the_pitchfork():
    rec = rc.record(PITCHFORK, FAST)
    assert rec["schema"] == rc.SCHEMA and rec["spec_sha256"] and rec["implementation_sha256"]
    t = rec["targets"]
    assert t["symmetric-write"]["status"].startswith("reached") and t["threshold-write"]["status"].startswith("reached")
    assert t["phase-locking"]["status"] == "obstructed" and "law_constant" not in t["symmetric-write"]
    comp = rec["completion"]
    assert comp["status"] == "derived" and len(comp["write_points"]) == 1
    law = comp["write_points"][0]["law"]
    assert law["status"] == "derived" and law["write_class"] == "linear-stage write"
    assert law["hold"]["switching_fraction"] < 0.05           # the written states retain the choice
    assert rec["retention"][0]["write_class"] == "linear-stage write"
    assert rec["attractor"] == {"status": "stable states", "stable_states": 2, "seconds": rec["attractor"]["seconds"]}
    assert rec["summary"].startswith("reaches symmetric-write, threshold-write")
    md = rc.markdown([rec])
    assert "| --- |" in md and "linear-stage write" in md and "2 stable states" in md


def test_the_record_command_isolates_its_steps_and_reports_a_stopped_one(tmp_path):
    specs = tmp_path / "specs"
    specs.mkdir()
    (specs / "pitchfork.json").write_text(json.dumps(PITCHFORK), encoding="utf-8")
    (specs / "survey.json").write_text(json.dumps({"carriers": []}), encoding="utf-8")   # not a specification: skipped
    out = tmp_path / "rec"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "record", str(specs), "--out-dir", str(out), "--no-law",
           "--budget", "240", "--completion-budget", "1", "--attractor-budget", "240", "--samples", "400",
           "--trajectories-completion", "300"]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    assert res.returncode == 0, res.stderr
    summary = json.loads(res.stdout)
    assert summary["records"] == 1
    rows = [json.loads(l) for l in (out / "records.jsonl").read_text(encoding="utf-8").splitlines() if l.strip()]
    rec = rows[0]
    assert rec["schema"] == rc.SCHEMA and rec["targets"]["symmetric-write"]["status"].startswith("reached")
    assert rec["write_points"]["status"] == "derived"
    # the completion child was stopped at one second; the record keeps the other steps and says why
    assert rec["completion"]["status"] == "not derived" and "time limit" in rec["completion"]["note"]
    assert rec["attractor"]["status"] == "stable states"
    md = (out / "records.md").read_text(encoding="utf-8")
    assert "| Specification |" in md and "not derived" in md


def test_the_attractor_step_runs_the_blow_up_law_for_a_slow_escape():
    growth = json.loads((ROOT / "examples" / "dynamics" / "hyperbolic_growth.json").read_text(encoding="utf-8"))
    out = rc.step(growth, "attractor")
    assert out["status"] == "measured" and out["stable_states"] == 0
    assert out["motion"]["kind"] in ("drift", "unbounded")
    assert out["laws"]["kind"] in ("blow-up", "growth"), out["laws"]


def test_the_table_of_the_write_targets_is_a_view_of_the_records(tmp_path):
    records = rc.load_records(ROOT / "docs" / "records" / "quantum_semiclassical.jsonl")
    assert len(records) == 31 and all(r["schema"] == rc.SCHEMA for r in records)
    rows = rc.predictions_rows(records)
    assert [r["spec"] for r in rows] == [Path(r["source"]).name for r in records]
    assert rows[0]["targets"] is records[0]["targets"] and rows[0]["rule"] == "write point"
    counts = rc.predictions_counts(records)
    assert counts["specifications"] == 31 and counts["reached"]["symmetric-write"] + counts["obstructed_symmetric"] <= 31
    assert set(counts["write_classes"]) <= {"linear-stage write", "equilibrium write", "nonlinear stage"}
    assert "kerr_cat_kpo" in counts["write_classes"]["equilibrium write"]
    md = rc.predictions_markdown(records, intro="Intro sentence.")
    assert md.startswith("# The write rules on the recorded specifications\n\nIntro sentence. Derived from the records")
    assert md.count("\n| ") == 31 + 2 and "equilibrium write: 51 % change state" in md
    # the command writes the same table from a records file, with no derivation
    out = tmp_path / "table"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "record", "--from-records",
           str(ROOT / "docs" / "records" / "quantum_semiclassical.jsonl"), "--out-dir", str(out)]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    assert res.returncode == 0, res.stderr
    assert json.loads(res.stdout)["specifications"] == 31
    assert (out / "quantum_semiclassical.md").read_text(encoding="utf-8").count("\n| ") == 33
    table = json.loads((out / "quantum_semiclassical.json").read_text(encoding="utf-8"))
    assert table["derived_from"].endswith("quantum_semiclassical.jsonl") and len(table["rows"]) == 31
    assert table["counts"] == counts


def test_the_cells_of_the_table_name_the_status_the_law_and_the_class():
    assert rc.target_cell({"status": "reached", "law_constant": {"constant": 1.3, "stderr": 0.1, "expected": 1.331}}) == \
        "reached; law constant 1.300 ± 0.100 (expected 1.331)"
    assert rc.target_cell({"status": "not derived", "note": "stopped at the time limit of 300 s"}) == \
        "not derived (stopped at the time limit of 300 s)"
    assert rc.target_cell({"status": "obstructed", "obstruction_short": "fold; no cusp"}) == "obstructed: fold; no cusp"
    law = {"status": "derived", "certified": False, "write_class": "equilibrium write",
           "hold": {"hold_time": 56.6, "switching_fraction": 0.51}}
    assert rc.completion_cell({"write_points": [{"completion": "leading order", "law": law}]}) == \
        "leading order: derived, not certified (equilibrium write: 51 % change state during a hold of 56.6)"
    assert rc.completion_cell({"status": "not derived", "note": "stopped at the time limit of 1 s"}) == \
        "not run (stopped at the time limit of 1 s)"
    assert rc.completion_cell({"status": "derived", "write_points": []}) == "no write point along the control"


def test_a_run_restricted_to_some_steps_records_the_others_as_not_run(tmp_path):
    specs = tmp_path / "specs"
    specs.mkdir()
    (specs / "pitchfork.json").write_text(json.dumps(PITCHFORK), encoding="utf-8")
    out = tmp_path / "rec"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "record", str(specs), "--out-dir", str(out), "--no-law",
           "--steps", "write_points,attractor", "--no-isolate"]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    assert res.returncode == 0, res.stderr
    rec = rc.load_records(out / "records.jsonl")[0]
    assert rec["write_points"]["status"] == "derived" and rec["attractor"]["status"] == "stable states"
    assert rec["completion"] == {"status": "not run"} and rec["targets"]["symmetric-write"] == {"status": "not run"}
    assert (out / "predictions.md").exists() and json.loads((out / "predictions.json").read_text(encoding="utf-8"))["counts"]["specifications"] == 1
