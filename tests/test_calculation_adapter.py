from copy import deepcopy
from dataclasses import replace
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

pytest.importorskip("sympy")
from fieldbridge.calculation_adapter import attach_calculation, canonical_source_equation
from fieldbridge.constructor import construct_transfer

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples/calculated_transfer"
DATA = EXAMPLES / "data"
QUERY = "radial Brownian Ito diffusion stochastic generator"


def request(name="squared_signal"):
    return json.loads((EXAMPLES / f"{name}.json").read_text())


def construct(req=None, data=DATA, query=QUERY, target="stochastic_dynamics"):
    return construct_transfer(query, target, data_dir=data, include_hyperion=False,
                              calculation_request=request() if req is None else req)


def test_constructor_retrieves_emits_and_calculates():
    result = construct()
    assert result.calculation["binding"]["record_id"] == "radial_ito_source"
    assert result.calculation["status"] == "calculated"
    assert result.calculation["report"]["target"]["drift"] == "-2*alpha*y + 2*theta + 2"
    assert result.calculation["report"]["candidate"]["residual_d_phi"] == "1"
    assert result.validation_gates["source_record_linked"]
    assert result.validation_gates["local_generator_identity"]
    assert all(not result.validation_gates[key] for key in (
        "field_blind_retrieval", "independent_target_validation", "prospective_prediction"))
    assert "Mean response:" in result.predictions[-1]


def test_emitted_source_comes_from_record_not_detached_example(tmp_path):
    data = tmp_path / "data"
    shutil.copytree(DATA, data)
    path = data / "index/core_examples.json"
    records = json.loads(path.read_text())
    records[0]["calculation_source"]["drift"] = "(theta+1/2)/x - 3*alpha*x"
    records[0]["equations"][0] = canonical_source_equation(records[0]["calculation_source"])
    path.write_text(json.dumps(records))
    result = construct(data=data)
    assert result.calculation["report"]["target"]["drift"] == "-6*alpha*y + 2*theta + 2"
    assert result.calculation["binding"]["source_record_sha256"] != construct().calculation["binding"]["source_record_sha256"]


def test_changing_map_changes_the_derived_result():
    req = request()
    req.update(state_map="2*x", inverse_map="y/2")
    req.pop("candidate")
    result = construct(req)
    assert result.calculation["status"] == "calculated"
    assert result.calculation["report"]["target"]["variance"] == "4"
    assert not result.validation_gates["omission_control_detected"]


def test_quantum_hamiltonian_is_retrieved_and_observable_is_supplied():
    result = construct(request("quantum_signal"), query="quantum Hamiltonian spin commutator", target="quantum_dynamics")
    assert result.calculation["status"] == "calculated"
    assert result.calculation["report"]["observable_dimension"] == 2
    assert result.validation_gates["finite_quantum_closure"]
    assert not result.validation_gates["local_generator_identity"]


@pytest.mark.parametrize("key,value", [
    ("source_record_id", "not_retrieved"), ("target_field", "different"),
    ("kind", "unsupported"), ("schema", "future/99"), ("assumptions", []),
    ("drift", "0"), ("state_map", "__import__('os')"), ("inverse_map", "y"),
    ("question", ""), ("observable", [[1]])])
def test_invalid_or_unsupported_correspondence_refuses_instead_of_guessing(key, value):
    req = request()
    req[key] = value
    result = construct(req)
    assert result.calculation["status"] == "refused"
    assert not result.predictions
    assert not result.validation_gates["local_generator_identity"]


def test_false_annotation_cannot_be_attached_to_different_stored_equation():
    baseline = construct_transfer(QUERY, "stochastic_dynamics", data_dir=DATA, include_hyperion=False)
    match = baseline.translation.matches[0]
    source = deepcopy(match.record.calculation_source)
    source["noise"] = "2"
    changed = replace(match, record=replace(match.record, calculation_source=source))
    baseline = replace(baseline, translation=replace(baseline.translation, matches=[changed]))
    result = attach_calculation(baseline, request())
    assert result.calculation["status"] == "refused"
    assert "does not match" in result.calculation["reason"]


def test_duplicate_retrieved_ids_refuse_ambiguous_binding():
    baseline = construct_transfer(QUERY, "stochastic_dynamics", data_dir=DATA, include_hyperion=False)
    baseline = replace(baseline, translation=replace(baseline.translation, matches=baseline.translation.matches*2))
    assert attach_calculation(baseline, request()).calculation["status"] == "refused"


def test_ordinary_retrieved_equation_does_not_gain_a_fabricated_model():
    baseline = construct_transfer("Brownian probability flow", "stochastic_optimization",
                                  data_dir=ROOT / "data", include_hyperion=False)
    req = request()
    req.update(source_record_id=baseline.translation.matches[0].record.record_id, target_field="stochastic_optimization")
    result = attach_calculation(baseline, req)
    assert result.calculation["status"] == "refused"
    assert "lacks a typed source model" in result.calculation["reason"]


def test_cli_writes_spec_and_clears_stale_success_on_refusal(tmp_path):
    correspondence = tmp_path / "correspondence.json"
    correspondence.write_text(json.dumps(request()))
    output = tmp_path / "result"
    command = [sys.executable, "-B", "-m", "fieldbridge", "--data-dir", str(DATA),
               "construct", QUERY, "--to", "stochastic_dynamics", "--no-hyperion", "--calculate",
               "--correspondence", str(correspondence), "--out-dir", str(output)]
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert "Derived Consequences" in result.stdout
    assert "Preserved Contract" not in result.stdout
    assert json.loads((output / "construction_spec.json").read_text())["drift"] == "(theta + 1/2)/x - alpha*x"
    bad = request()
    bad["source_record_id"] = "absent"
    correspondence.write_text(json.dumps(bad))
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    assert result.returncode == 2
    assert json.loads((output / "run_status.json").read_text())["status"] == "refused"
    assert json.loads((output / "construction_spec.json").read_text()) is None
    assert json.loads((output / "calculation.json").read_text())["status"] == "refused"


def test_without_calculate_keeps_the_original_proposal_path():
    result = construct_transfer(QUERY, "stochastic_dynamics", data_dir=DATA, include_hyperion=False)
    assert result.calculation is None
    assert len(result.predictions) == 3


def test_cli_requires_explicit_calculation_inputs(tmp_path):
    result = subprocess.run([sys.executable, "-B", "-m", "fieldbridge", "construct", QUERY,
                             "--to", "stochastic_dynamics", "--calculate"],
                            cwd=ROOT, text=True, capture_output=True)
    assert result.returncode != 0
    assert "requires --correspondence and --out-dir" in result.stderr


def test_no_model_can_be_injected_by_correspondence():
    req = request()
    req["parameters"] = {"theta": "positive"}
    result = construct(req)
    assert result.calculation["status"] == "refused"
    assert "unsupported" in result.calculation["reason"]
