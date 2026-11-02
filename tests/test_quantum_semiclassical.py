"""The quantum-computing elements written as semiclassical drift equations (examples/quantum/semiclassical), and the
write rules over a directory of drift specifications (memory predictions)."""
from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

from fieldbridge.memory import spec as mspec

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "examples" / "quantum" / "semiclassical"
CARRIERS = [
    "bright_state_resonator", "cim_dopo", "cim_pair", "cim_triangle", "dissipative_cat", "flux_qubit_tilt", "fluxonium",
    "ion_duffing", "ion_phonon_laser", "jba_resonator", "jpc_nondegenerate", "kerr_cat_kpo", "kpo_pair", "kpo_triangle",
    "magnon_kerr", "nuclear_bath", "optomech_bistable", "optomech_selfosc", "parametric_amplifier",
    "phase_qubit_rfsquid", "po_amplitude_oscillator", "pplo_oscillator", "qfp_latch", "quantum_vdp", "rfsquid_annealer",
    "rfsquid_pair", "rfsquid_triangle", "rydberg_bistability", "snail_po", "spin_ensemble_cavity", "transmon_driven",
]


def test_the_folder_holds_the_thirty_one_elements_of_the_census():
    specs = [p for p in FOLDER.glob("*.json") if "schema" in json.loads(p.read_text(encoding="utf-8"))]
    assert sorted(p.stem for p in specs) == CARRIERS   # census.json is the survey, not a specification
    readme = (FOLDER / "README.md").read_text(encoding="utf-8")
    assert all(f"`{name}.json`" in readme for name in CARRIERS)


@pytest.mark.parametrize("name", CARRIERS)
def test_every_element_loads_and_states_its_source_and_assumptions(name):
    path = FOLDER / f"{name}.json"
    spec = json.loads(path.read_text(encoding="utf-8"))
    assert spec["schema"] == "fieldbridge-memory/1" and spec["kind"] == "equations"
    assert re.search(r"(19|20)\d\d", spec["provenance"]["source"]), "a dated published source"
    assert spec["assumptions"] and spec["field"] and spec["closure"] and spec["observable"]
    assert spec["question"]
    real = mspec.load(path)
    assert real.control == spec["control"]["name"] and real.control in spec["parameters"]
    lo, hi = spec["control"]["range"]
    assert lo < hi
    assert set(spec.get("material", [])) <= set(spec["parameters"]) and real.control not in spec.get("material", [])
    assert real.carrier.dim == len(spec["carrier"]["variables"]) >= 1



def test_predictions_derive_the_writes_of_drift_specifications_and_complete_an_obstruction(tmp_path):
    """Two drift specifications: the pitchfork, whose symmetric write is reached, and a subcritical pitchfork saturated
    at fifth order, whose symmetric write is obstructed and completed on the distant branch (the choice is made in the
    linear stage; the positive cubic only decides where the state jumps)."""
    specs = tmp_path / "specs"
    specs.mkdir()
    pitch = json.loads((ROOT / "examples" / "memory" / "pitchfork.json").read_text(encoding="utf-8"))
    (specs / "pitchfork.json").write_text(json.dumps(pitch), encoding="utf-8")
    sub = dict(pitch)
    sub.pop("potential", None)
    var = sub["carrier"]["variables"][0]
    sub["name"] = "subcritical pitchfork saturated at fifth order"
    sub["drift"] = {var: f"{sub['control']['name']}*{var} + {var}**3 - {var}**5"}
    sub["control"] = {"name": sub["control"]["name"], "range": [-0.5, 0.5]}
    (specs / "subcritical.json").write_text(json.dumps(sub), encoding="utf-8")
    out = tmp_path / "pred"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "predictions", str(specs), "--out-dir", str(out),
           "--no-law", "--budget", "300", "--completion-budget", "600", "--samples", "300",
           "--trajectories-completion", "200"]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    assert res.returncode == 0, res.stderr
    summary = json.loads(res.stdout)
    assert summary["drift"] == 2 and summary["specifications"] == 2
    report = json.loads((out / "predictions.json").read_text(encoding="utf-8"))
    rows = {r["spec"]: r for r in report["rows"]}
    assert set(rows) == {"pitchfork.json", "subcritical.json"} and all(r["rule"] == "write point" for r in rows.values())
    p = rows["pitchfork.json"]["targets"]
    assert p["symmetric-write"]["status"].startswith("reached") and p["threshold-write"]["status"].startswith("reached")
    assert p["phase-locking"]["status"] == "obstructed" and "law_constant" not in p["symmetric-write"]  # --no-law
    assert "completion" not in rows["pitchfork.json"]  # nothing to complete when the write is reached
    s = rows["subcritical.json"]
    assert s["targets"]["symmetric-write"]["status"] == "obstructed"
    assert s["targets"]["threshold-write"]["status"].startswith("reached")
    points = s["completion"]["write_points"]
    kinds = {e["completion"] for e in points}
    assert "distant branch" in kinds, points
    branch = next(e for e in points if e["completion"] == "distant branch")
    assert branch["law"]["status"] == "derived" and isinstance(branch["law"]["certified"], bool)
    assert abs(branch["write_point"][sub["control"]["name"]]) < 0.05  # the pitchfork of the subcritical form sits at zero
    md = (out / "predictions.md").read_text(encoding="utf-8")
    assert "| Specification | Field | Symmetric write | Threshold write | Phase locking | Completion at the obstruction |" in md
    assert "distant branch: derived" in md and "fold" in md


def test_a_child_that_outruns_its_time_limit_keeps_the_steps_it_had_written(tmp_path):
    """Every step runs in a child of its own under its own wall-clock limit: thirty seconds for each write target,
    which take a few, and one second for a completion of forty thousand trajectories, which takes minutes. The
    completion child is stopped; the row keeps the three targets, says why the completion is missing, and the run
    goes on."""
    specs = tmp_path / "specs"
    specs.mkdir()
    pitch = json.loads((ROOT / "examples" / "memory" / "pitchfork.json").read_text(encoding="utf-8"))
    pitch.pop("potential", None)
    var = pitch["carrier"]["variables"][0]
    pitch["name"] = "subcritical pitchfork saturated at fifth order"
    pitch["drift"] = {var: f"{pitch['control']['name']}*{var} + {var}**3 - {var}**5"}
    pitch["control"] = {"name": pitch["control"]["name"], "range": [-0.5, 0.5]}
    (specs / "subcritical.json").write_text(json.dumps(pitch), encoding="utf-8")
    out = tmp_path / "pred"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "predictions", str(specs), "--out-dir", str(out),
           "--no-law", "--budget", "30", "--completion-budget", "1", "--samples", "40000",
           "--trajectories-completion", "40000"]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    assert res.returncode == 0, res.stderr
    row = json.loads((out / "predictions.json").read_text(encoding="utf-8"))["rows"][0]
    assert row["targets"]["symmetric-write"]["status"] == "obstructed"
    assert row["targets"]["threshold-write"]["status"].startswith("reached")
    assert "stopped at the time limit" in row["completion"]["note"] and row["completion"]["write_points"] == []
    assert "stopped at the time limit" in (out / "predictions.md").read_text(encoding="utf-8")


# ------------------------------------------------------------------------------------------------ the survey and Module 15
def _survey_script(name: str):
    spec = importlib.util.spec_from_file_location(f"survey_{name}", FOLDER / "survey" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _survey_results():
    return {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in sorted((FOLDER / "survey" / "results").glob("*.json"))}


def _records():
    lines = (ROOT / "docs" / "records" / "quantum_semiclassical.jsonl").read_text(encoding="utf-8").splitlines()
    return {Path(r["source"]).stem: r for r in map(json.loads, lines)}


def test_the_records_and_the_survey_are_of_the_current_specifications_and_of_one_code():
    records, results = _records(), _survey_results()
    assert sorted(records) == CARRIERS and sorted(results) == CARRIERS
    for name in CARRIERS:
        spec = json.loads((FOLDER / f"{name}.json").read_text(encoding="utf-8"))
        assert records[name]["spec_sha256"] == mspec.input_hash(spec), name
    assert len({r["implementation_sha256"] for r in records.values()}) == 1
    assert len({r["implementation_sha256"] for r in results.values()}) == 1


def test_module_15_states_the_counts_of_the_survey():
    report = _survey_script("report")
    results = _survey_results()
    text = " ".join((ROOT / "docs" / "tutorial" / "35_memory_quantum_hardware.md").read_text(encoding="utf-8").split())
    cells = {t: [report.derivation_cell(t, r["derivations"][t]) for r in results.values()] for t in report.EXPECTED}

    def known(t):
        return sum(c.startswith("known law") for c in cells[t])

    def reached(t):
        return sum(c.startswith(("known law", "law differs", "reached")) for c in cells[t])

    assert f"reached in {reached('symmetric-write')} of the 31 models. In {known('symmetric-write')} of them" in text
    assert f"reached in {reached('threshold-write')} models, and in {known('threshold-write')} of them" in text
    laws = [x for r in results.values() for x in report.law_rows(r["completion"]) if x["status"].startswith("derived")]
    certified = sum(bool(x.get("certified")) for x in laws)
    assert f"{certified} laws at obstructions are certified and {len(laws) - certified} are derived" in text
    changed = sum(any(s["reaches"] for s in r["changes"]["singles"]) for r in results.values())
    assert f"In {changed} of the 31 models a single change" in text


@pytest.mark.parametrize("name", ["phase_qubit_rfsquid", "magnon_kerr"])
def test_the_write_points_counted_once_at_a_fold_are_folds_in_the_records_and_the_survey(name):
    """At a fold the growth along the branch diverges; it once took every eigenvalue of these two elements into the
    critical count, and the records read a degenerate critical mode."""
    record = _records()[name]
    assert record["write_points"]["points"]
    for point in record["write_points"]["points"]:
        assert point["kind"].startswith("saddle-node") and point["critical_dimension"] == 1
    assert {e["completion"] for e in record["completion"]["write_points"]} == {"second parameter"}
    assert {e["completion"] for e in _survey_results()[name]["completion"]["write_points"]} == {"second parameter"}


def test_the_phase_qubit_reaches_the_symmetric_write_through_its_cusp_and_the_threshold_write_along_the_flux():
    targets = _records()["phase_qubit_rfsquid"]["targets"]
    assert targets["symmetric-write"]["status"] == "reached after unfolding" and targets["symmetric-write"]["word"] == "SCRURKL"
    assert targets["threshold-write"]["status"] == "reached" and targets["threshold-write"]["word"] == "CRKL"


def test_the_start_of_the_check_lowers_the_constant_of_a_weakly_damped_element():
    """The linear stage alone, integrated exactly: from rest 4/sqrt(r) before the crossing, as the check starts, the
    ringing of the phase particle lowers the constant where the damping is weak; a start twice as early removes it."""
    start = _survey_script("phase_qubit_start")
    assert start.linear_constant(0.1) == pytest.approx(1.169, abs=0.005)
    assert start.linear_constant(0.2) == pytest.approx(1.311, abs=0.005)
    assert start.linear_constant(0.1, before=8.0) == pytest.approx(1.332, abs=0.005)
    assert start.linear_constant(0.2, before=8.0) == pytest.approx(1.324, abs=0.005)
    assert start.linear_constant(1.0) == pytest.approx(start.linear_constant(1.0, before=8.0), abs=0.001)
    # with no noise on the phase itself the law of one coordinate holds to 0.2 %
    assert start.linear_constant(0.2, before=8.0, noise_on_phase=False) == pytest.approx(start.math.pi ** 0.25, rel=0.002)
