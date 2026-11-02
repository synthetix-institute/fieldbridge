import copy
import json
from pathlib import Path

import pytest

from fieldbridge.memory import learning as ln
from fieldbridge.memory import repair_families as rf


def task(name="hidden_exchange"):
    return json.loads((Path(__file__).resolve().parents[1] / "examples/learning" / (name + ".json")).read_text())


def trained(tmp_path, name="hidden_exchange"):
    library = rf.FamilyLibrary(tmp_path / "families.json")
    t = task(name)
    r = ln.learn(t)
    assert r["status"] == "verified_extension"
    assert library.abstract(t, r, "training")
    return library


def test_mixed_terms_require_an_added_coordinate():
    t = task()
    t["model"]["drift"]["x"] = "-a*x+b*x*y"
    tower = rf.omitted_tower(t, 2)
    assert tower["closed"]
    assert tower["added"] == [[0, 1]]
    result = rf.construct(t, strategy="symbolic")
    assert result["status"] == "verified_extension", result
    assert result["candidate"]["coordinates"]["retained_0"] == "y"


def test_budget_exhaustion_is_not_closure():
    t = task()
    t["model"]["drift"]["x"] = "-a*x+b*x*y"
    assert not rf.omitted_tower(t, 0)["closed"]
    result = rf.construct(t, strategy="symbolic", policy=ln.Policy(max_added_coordinates=0))
    assert result["status"] != "verified_extension"


def test_revised_hint_is_skipped_and_fresh_policy_remains_eligible(tmp_path):
    library = trained(tmp_path)
    t = task()
    policy = ln.Policy(max_candidates=1)
    before, proposals = rf._proposals(t, library, policy)
    assert proposals
    library.revise("failed-hint", t, "omitted_state_variable", "budget-limited hint", policy=policy)
    after, proposals = rf._proposals(t, library, policy)
    assert not proposals and "guidance_skipped" in after
    assert rf._proposals(t, library, ln.Policy())[1]
    t["withheld"][0]["id"] += "-renamed"
    assert not rf._proposals(t, library, policy)[1]


def test_frozen_library_exposes_only_copies(tmp_path):
    library = trained(tmp_path)
    frozen = rf.FamilyLibrary(library.path, frozen=True)
    before = frozen.digest()
    frozen.data["families"].clear()
    frozen.families().clear()
    assert frozen.digest() == before
    with pytest.raises(RuntimeError):
        frozen.revise("x", task(), "omitted_state_variable", "test")


def test_family_proposals_obey_degree_and_term_limits(tmp_path):
    library = trained(tmp_path, "quintic_saturation")
    policy = ln.Policy(max_degree=1, max_terms=1)
    report = rf.construct(task("quintic_saturation"), library, policy)
    assert report["status"] != "verified_extension"
    assert report.get("rejected_proposals")


def test_family_runner_consumes_withheld_tests(tmp_path):
    library = trained(tmp_path)
    evaluation = ln.MoveLibrary(tmp_path / "tests.sqlite")
    assert rf.construct(task(), library, evaluation_library=evaluation)["status"] == "verified_extension"
    again = task()
    again["candidate_drift"]["u"] = "0"
    report = rf.construct(again, library, evaluation_library=evaluation)
    assert report["status"] == "withheld_already_used"
    assert "withheld" not in report


def test_no_family_can_be_seeded_with_an_unverified_claim(tmp_path):
    library = rf.FamilyLibrary(tmp_path / "families.json")
    assert library.abstract(task(), {"status": "verified_extension"}, "fake") is None
    assert not library.families()


def test_departure_to_repair_to_reuse_is_connected(tmp_path):
    t = task()
    evaluations = ln.MoveLibrary(tmp_path / "evaluations.sqlite")
    library = rf.FamilyLibrary(tmp_path / "families.json")
    record = {"spec": t["model"], "targets": {"reduction": {
        "status": "obstructed", "learning_task": t}}}
    ln.retain_record(evaluations, record, source="source.json")
    rows = rf.repair_departures(evaluations, library)
    assert rows[0]["status"] == "verified_extension", rows
    assert rows[0]["retained_family"] == "omitted_state_variable"
    assert rf.repair_departures(evaluations, library)[0]["status"] == "already_processed"
    other = copy.deepcopy(t)
    other["model"]["drift"]["x"] = "-a*x+2*b*y"
    fresh = rf.construct(other, rf.FamilyLibrary(library.path, frozen=True), evaluation_library=evaluations)
    assert fresh["status"] == "verified_extension", fresh
    assert fresh["selection_origin"] == "families"


def test_unstructured_departure_and_stochastic_task_remain_unresolved(tmp_path):
    evaluations = ln.MoveLibrary(tmp_path / "evaluations.sqlite")
    library = rf.FamilyLibrary(tmp_path / "families.json")
    ln.retain_record(evaluations, {"targets": {"write": {"status": "not derived"}}}, source="old.jsonl")
    assert rf.repair_departures(evaluations, library)[0]["status"] == "needs_equation_task"
    t = task()
    t["model"]["noise"] = .1
    ln.retain_record(evaluations, {"spec": t["model"], "targets": {"write": {
        "status": "obstructed", "learning_task": t}}}, source="stochastic.json")
    rows = rf.repair_departures(evaluations, library)
    assert any(r["status"] == "unsupported_or_invalid" for r in rows)
    assert not library.families()


def test_public_repair_command_reads_recorded_failures(tmp_path):
    from fieldbridge.cli import main
    t = task()
    registry = tmp_path / "registry.sqlite"
    evaluations = ln.MoveLibrary(registry)
    ln.retain_record(evaluations, {"spec": t["model"], "targets": {"reduction": {
        "status": "obstructed", "learning_task": t}}}, source="source.json")
    code = main(["memory", "repair", "--library", str(registry), "--families", str(tmp_path / "families.json"),
                 "--out-dir", str(tmp_path / "out")])
    assert code == 0
    assert json.loads((tmp_path / "out/repair.json").read_text())["campaigns"][0]["status"] == "verified_extension"


def test_family_and_direct_solver_share_numerical_failure_handling(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise ArithmeticError("solver did not converge")
    monkeypatch.setattr(ln, "_calculate_cases", fail)
    result = rf.construct(task(), strategy="symbolic")
    assert result["status"] == "numerical_failure"
    assert "withheld" not in result


def test_source_placeholder_names_do_not_collide():
    import sympy as sp
    x = sp.Symbol("_z0", real=True)
    rem, _, tags = rf.reduce_in_subalgebra(x**2, [x], [x], [])
    assert not rem.free_symbols.intersection({x})
    assert rem == tags[0]**2


def test_interrupted_campaign_is_not_reported_as_a_completed_repair(tmp_path):
    from fieldbridge.cli import main
    t = task()
    registry = tmp_path / "registry.sqlite"
    evaluations = ln.MoveLibrary(registry)
    record = {"spec": t["model"], "targets": {"reduction": {"status": "obstructed", "learning_task": t}}}
    key = ln.retain_record(evaluations, record, source="source.json")[0]["id"]
    assert evaluations.claim_campaign(key)
    assert not evaluations.claim_campaign(key)
    families = rf.FamilyLibrary(tmp_path / "families.json")
    assert rf.repair_departures(evaluations, families)[0]["status"] == "campaign_in_progress"
    assert main(["memory", "repair", "--library", str(registry), "--families", str(families.path),
                 "--out-dir", str(tmp_path / "out")]) == 2
