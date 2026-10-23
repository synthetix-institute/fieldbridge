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
