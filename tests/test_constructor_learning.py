import copy
import json
from pathlib import Path

import pytest
import sympy as sp

from fieldbridge.memory import learning as L

EX = Path(__file__).resolve().parents[1] / "examples/learning"


def task(name="quintic_saturation"):
    return json.loads((EX / f"{name}.json").read_text())


def test_derives_quintic_without_a_target_constant_and_reuses_it(tmp_path):
    library = L.MoveLibrary(tmp_path / "moves.sqlite")
    first = L.learn(task(), library)
    assert first["status"] == "verified_extension", first
    assert first["candidate"]["drift"]["u"] == "a*u - b*u**5"
    assert first["initial_residual"] == ["-b*x**5 + b*x**3"]
    assert first["withheld"]["max_error"] < 1e-7
    assert first["withheld"]["max_omission_error"] > .01
    second_task = task()
    second_task["model"]["parameters"]["k"] = .2
    second_task["model"]["drift"]["x"] = "a*x-k*x**3-2*b*x**5"
    second_task["candidate_drift"]["u"] = "a*u-k*u**3"
    second = L.learn(second_task, library)
    assert second["status"] == "verified_extension", second
    assert second["selection_origin"] == "learned"
    assert second["candidate_evaluations"] < first["candidate_evaluations"]
    assert second["candidate"]["drift"]["u"] == "a*u - 2*b*u**5 - k*u**3"
    assert second["candidate_sha256"] != first["candidate_sha256"]
    assert not first["novelty_established"]


def test_hidden_initial_state_is_retained_even_without_feedback(tmp_path):
    result = L.learn(task("hidden_exchange"), L.MoveLibrary(tmp_path / "moves.sqlite"))
    assert result["status"] == "verified_extension", result
    assert result["candidate"]["coordinates"] == {"u": "x", "retained_0": "y"}
    assert result["candidate"]["drift"] == {"u": "-a*u + b*retained_0", "retained_0": "-c*retained_0"}
    assert result["withheld"]["max_omission_error"] > .01


def test_missing_coordinate_stays_unresolved_when_attachments_forbidden():
    r = L.learn(task("hidden_exchange"), policy=L.Policy(max_added_coordinates=0))
    assert r["status"] == "unresolved"
    assert "candidate" not in r


def test_exact_existing_model_needs_no_new_move(tmp_path):
    t = task()
    t["candidate_drift"]["u"] = "a*u-b*u**5"
    library = L.MoveLibrary(tmp_path / "library.sqlite")
    r = L.learn(t, library)
    assert r["status"] == "already_closed"
    assert not library.recipes()


@pytest.mark.parametrize("edit", ["noise", "orthant", "periodic", "nonpolynomial", "duplicate", "unseen", "syntax"])
def test_unsupported_and_invalid_tasks_cannot_be_learned(tmp_path, edit):
    t = task()
    if edit == "noise":
        t["model"]["noise"] = .1
    elif edit == "orthant":
        t["model"]["carrier"]["kind"] = "orthant"
    elif edit == "periodic":
        t["model"]["carrier"]["kind"] = "torus"
    elif edit == "nonpolynomial":
        t["model"]["drift"]["x"] = "sin(x)"
    elif edit == "duplicate":
        t["withheld"][0] = copy.deepcopy(t["development"][0])
        t["withheld"][0]["id"] = "another-name-does-not-make-independent-data"
        t["withheld"][0]["duration"] = 4
    elif edit == "unseen":
        for row in t["withheld"]:
            row["parameters"] = t["development"][0]["parameters"]
    else:
        t["candidate_drift"]["u"] = "__import__('os').getcwd()"
    library = L.MoveLibrary(tmp_path / "library.sqlite")
    r = L.learn(t, library)
    assert r["status"] == "unsupported_or_invalid", r
    assert not library.recipes()


def test_independent_verifier_rejects_corrected_flag_with_wrong_equation():
    t = task()
    candidate = {"coordinates": t["coordinates"], "drift": t["candidate_drift"], "passed": True}
    assert not L.verify_identity(t, candidate)["passed"]


def test_cannot_change_observable_to_make_a_test_pass():
    with pytest.raises(L.ConstructionError, match="observable"):
        L.verify_identity(task(), {"coordinates": {"u": "x**2"}, "drift": {"u": "0"}})


def test_polynomial_map_closes_without_a_full_state_inverse():
    t = task()
    t["model"]["drift"]["x"] = "-a*x"
    t["coordinates"] = {"u": "x**2"}
    t["candidate_drift"] = {"u": "-a*u"}
    r = L.learn(t)
    assert r["status"] == "verified_extension", r
    assert r["candidate"]["drift"]["u"] == "-2*a*u"


def test_new_observable_is_derived_from_the_residual_not_copied_from_state():
    t = task("hidden_exchange")
    t["model"]["drift"]["x"] = "-a*x+b*y**2"
    r = L.learn(t)
    assert r["status"] == "verified_extension", r
    assert r["candidate"]["coordinates"] == {"u": "x", "retained_0": "y**2"}
    assert r["candidate"]["drift"]["retained_0"] == "-2*c*retained_0"


def test_no_learning_from_unexcited_preparations(tmp_path):
    t = task("hidden_exchange")
    for i, row in enumerate(t["development"]):
        row["initial"] = [i + .2, 0]
    library = L.MoveLibrary(tmp_path / "library.sqlite")
    r = L.learn(t, library)
    assert r["status"] == "insufficient_excitation", r
    assert not library.recipes()


def test_solver_failure_is_not_a_physical_obstruction(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise ArithmeticError("did not converge")
    monkeypatch.setattr(L, "_trajectory", fail)
    library = L.MoveLibrary(tmp_path / "library.sqlite")
    r = L.learn(task(), library)
    assert r["status"] == "numerical_failure"
    assert not library.recipes()


def test_withheld_failure_does_not_trigger_reselection_or_promotion(tmp_path, monkeypatch):
    actual = L._calculate_cases
    calls = []
    def calculate(task, candidate, cases, policy):
        calls.append([r["id"] for r in cases])
        result = actual(task, candidate, cases, policy)
        if cases[0]["id"].startswith("intervene"):
            result.update(passed=False, max_error=.5)
        return result
    monkeypatch.setattr(L, "_calculate_cases", calculate)
    library = L.MoveLibrary(tmp_path / "library.sqlite")
    r = L.learn(task(), library)
    assert r["status"] == "withheld_failed"
    assert len(calls) == 2
    assert not library.recipes()


def test_budget_exhaustion_does_not_become_absence():
    r = L.learn(task(), policy=L.Policy(max_candidates=1))
    assert r["status"] == "search_budget_exhausted"


def test_persistent_library_blocks_reusing_withheld_outcomes_for_a_new_trial(tmp_path):
    library = L.MoveLibrary(tmp_path / "moves.sqlite")
    assert L.learn(task(), library)["status"] == "verified_extension"
    t = task()
    t["candidate_drift"]["u"] = "a*u"
    for row in t["withheld"]:
        row["id"] += "-renamed"
        row["duration"] += 1
    r = L.learn(t, library)
    assert r["status"] == "withheld_already_used", r
    assert "withheld" not in r


def test_library_rejects_tampered_recipe(tmp_path):
    library = L.MoveLibrary(tmp_path / "moves.sqlite")
    assert L.learn(task(), library)["status"] == "verified_extension"
    with library.connect() as db:
        db.execute("UPDATE moves SET recipe=?", (json.dumps({"engine": L.ENGINE, "powers": [[99]]}),))
    assert not library.recipes()


def test_historical_failures_are_kept_without_promotion(tmp_path):
    library = L.MoveLibrary(tmp_path / "moves.sqlite")
    record = {"paper_id": "example", "targets": {
        "symmetric-write": {"status": "not derived", "obstruction": "SVD failed"},
        "threshold-write": {"status": "obstructed", "obstruction": "no fold"},
        "phase-locking": {"status": "reached", "law_constant": {"constant": 1.0}}}}
    rows = L.retain_record(library, record, source="historical.jsonl")
    assert len(rows) == 3
    assert rows[0]["kind"] == "numerical_or_execution_failure"
    assert not library.recipes()
    L.retain_record(library, record, source="historical.jsonl")
    with library.connect() as db:
        assert db.execute("SELECT count(*) FROM departures").fetchone()[0] == 3


def test_completion_departures_are_not_lost(tmp_path):
    library = L.MoveLibrary(tmp_path / "moves.sqlite")
    record = {"completion": {"status": "derived", "write_points": [
        {"law": {"status": "derived", "certified": False, "write_class": "nonlinear stage"}}]}}
    rows = L.retain_record(library, record, source="completion.jsonl")
    assert any(r["kind"] == "unexplained_law_departure_requires_convergence_checks" for r in rows)
    assert not library.recipes()


def test_historical_nan_is_preserved_explicitly_without_promoting(tmp_path):
    library = L.MoveLibrary(tmp_path / "moves.sqlite")
    record = {"targets": {"write": {"status": "reached", "law_constant": {"constant": float("nan")}}}}
    assert len(L.retain_record(library, record, source="legacy.jsonl")) == 1
    with library.connect() as db:
        stored = json.loads(db.execute("SELECT record FROM departures").fetchone()[0])
    assert stored["law"]["constant"] == {"nonfinite_value": "nan"}
    assert stored["nonfinite_fields_in_record"]
    assert not library.recipes()


def test_parameter_bounds_are_kept():
    t = task()
    t["withheld"][0]["parameters"]["b"] = -1
    assert L.learn(t)["status"] == "unsupported_or_invalid"


def test_nonfinite_input_is_retained_as_invalid(tmp_path):
    t = task()
    t["model"]["parameters"]["b"] = float("nan")
    library = L.MoveLibrary(tmp_path / "moves.sqlite")
    assert L.learn(t, library)["status"] == "unsupported_or_invalid"
    with library.connect() as db:
        assert db.execute("SELECT count(*) FROM attempts").fetchone()[0] == 1
    assert not library.recipes()


def test_numeric_spelling_does_not_make_a_new_intervention():
    t = task()
    row = copy.deepcopy(t["development"][0])
    row["id"] = "looks-different"
    row["parameters"]["b"] = 1.0
    t["withheld"][0] = row
    assert L.learn(t)["status"] == "unsupported_or_invalid"


def test_record_cli_collects_outcomes_in_the_learning_library(tmp_path, monkeypatch):
    from fieldbridge.cli import main
    from fieldbridge.memory import cli, record
    source = tmp_path / "source.json"
    model = task()["model"]
    source.write_text(json.dumps(model))
    monkeypatch.setattr(cli, "_record_one", lambda *a, **k: record.assemble(
        model, {"symmetric-write": {"status": "not derived", "note": "time limit"}}, source=str(source)))
    db_path = tmp_path / "moves.sqlite"
    assert main(["memory", "record", str(source), "--no-isolate", "--out-dir", str(tmp_path / "out"),
                 "--learning-library", str(db_path)]) == 0
    library = L.MoveLibrary(db_path)
    with library.connect() as db:
        rows = [json.loads(r[0]) for r in db.execute("SELECT record FROM departures")]
    assert any(r["kind"] == "numerical_or_execution_failure" for r in rows)
    assert all(r["model"] == model for r in rows)
    assert not library.recipes()


def test_parameter_singularity_is_not_hidden_by_numeric_test_values():
    t = task()
    t["model"]["drift"]["x"] = "x/a"
    assert L.learn(t)["status"] == "unsupported_or_invalid"


def test_cli_is_integrated(tmp_path):
    from fieldbridge.cli import main
    code = main(["memory", "learn", str(EX / "quintic_saturation.json"),
                 "--library", str(tmp_path / "moves.sqlite"), "--out-dir", str(tmp_path / "report")])
    assert code == 0
    assert json.loads((tmp_path / "report/learning.json").read_text())["status"] == "verified_extension"
