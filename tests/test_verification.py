import copy
import json
from pathlib import Path
import subprocess
import sys

import pytest

sp = pytest.importorskip("sympy")
from fieldbridge.verification import ConstructionError, expression, verify_construction, run_file


ROOT = Path(__file__).resolve().parents[1]


def example(name):
    return json.loads((ROOT / "examples" / "construction" / f"{name}.json").read_text())


def test_square_map_derives_missing_drift_and_rejects_naive_candidate():
    report = verify_construction(example("ito_square"))
    assert report["status"] == "verified_local_generator_identity"
    assert report["residual_coefficients"] == ["0", "0"]
    assert report["candidate"]["passes_local_identity"] is False
    assert report["candidate"]["residual_d_phi"] == "1"
    assert report["candidate"]["residual_d2_phi"] == "0"
    assert report["omission_control"]["residual_d_phi"] == "1"
    assert report["source_alignment_verified"] is False
    assert report["novelty_established"] is False


def test_affine_map_needs_no_extra_drift():
    report = verify_construction(example("affine_transfer"))
    assert report["candidate"]["passes_local_identity"]
    assert not report["omission_control"]["detects_omission"]


def test_second_order_coefficient_is_checked_even_if_mean_is_right():
    spec = example("ito_square")
    spec["candidate"] = {"drift": "2*theta+2-2*alpha*y", "variance": "0"}
    result = verify_construction(spec)["candidate"]
    assert not result["passes_local_identity"]
    assert result["residual_d_phi"] == "0"
    assert result["residual_d2_phi"] != "0"


def test_mean_prediction_satisfies_derived_moment_equation():
    report = verify_construction(example("ito_square"))
    names = {n: sp.Symbol(n, positive=True) for n in ("theta", "alpha")}
    names.update({n: sp.Symbol(n, real=True) for n in ("m0", "t")})
    mean = expression(report["mean_prediction"]["expression"], names)
    assert sp.simplify(sp.diff(mean, names["t"]) - (2*names["theta"]+2-2*names["alpha"]*mean)) == 0
    assert sp.simplify(mean.subs(names["t"], 0) - names["m0"]) == 0


@pytest.mark.parametrize("value", ["0=1", "1=2", "__import__('os')", "x.real", "unknown+x", "[x]", "1/0"])
def test_expression_rejects_contradictions_and_unsupported_syntax(value):
    with pytest.raises(ConstructionError):
        expression(value, {"x": sp.Symbol("x")})


def test_inverse_must_be_valid_on_declared_domain():
    spec = example("ito_square")
    spec["domain"] = "real"
    with pytest.raises(ConstructionError, match="inverse"):
        verify_construction(spec)


def test_candidate_cannot_leak_source_coordinate():
    spec = example("ito_square")
    spec["candidate"]["drift"] = "x"
    with pytest.raises(ConstructionError, match="target"):
        verify_construction(spec)


def test_quantum_algorithm_finds_correlation_and_exact_closure():
    report = verify_construction(example("quantum_correlations"))
    assert report["observable_dimension"] == 2
    assert all(report["closed_identities"])
    assert report["evolution_matrix"] == [["0", "1"], ["-4*g**2", "0"]]
    assert report["omission_control"]["single_observable_closed"] is False


def test_closure_predicts_full_unitary_evolution_without_fitting():
    spec = example("quantum_correlations")
    spec["parameters"] = {}
    spec["hamiltonian"] = [[str(x).replace("g", "1") for x in row] for row in spec["hamiltonian"]]
    report = verify_construction(spec)
    basis = [sp.Matrix([[expression(v, {"I": sp.I}) for v in row] for row in A]) for A in report["basis"]]
    G = sp.Matrix([[expression(v, {}) for v in row] for row in report["evolution_matrix"]])
    H = sp.diag(1, -1, -1, 1)
    # |+y> on the first spin and |0> on the second has a nonzero hidden correlation.
    psi = sp.Matrix([1, 0, sp.I, 0]) / sp.sqrt(2)
    rho = psi * psi.conjugate().T
    t = sp.pi / 8
    U = (-sp.I * H * t).exp()
    full = sp.trace(U * rho * U.conjugate().T * basis[0])
    initial = sp.Matrix([sp.trace(rho * A) for A in basis])
    reduced = ((G * t).exp() * initial)[0]
    assert sp.simplify(sp.expand_complex(full - reduced)) == 0
    assert sp.simplify(full - initial[0]) != 0


def test_zero_coupling_restores_one_observable_closure():
    spec = example("quantum_correlations")
    spec["hamiltonian"] = [[0]*4 for _ in range(4)]
    report = verify_construction(spec)
    assert report["observable_dimension"] == 1
    assert report["omission_control"]["single_observable_closed"]


def test_nonhermitian_hamiltonian_is_rejected():
    spec = example("quantum_correlations")
    spec["hamiltonian"][0][1] = 1
    with pytest.raises(ConstructionError, match="Hermitian"):
        verify_construction(spec)


def test_output_hash_tracks_the_actual_model(tmp_path):
    spec = example("ito_square")
    source = tmp_path / "input.json"
    source.write_text(json.dumps(spec))
    report = run_file(source, tmp_path / "out")
    assert json.loads((tmp_path / "out" / "calculation.json").read_text()) == report
    changed = copy.deepcopy(spec)
    changed["drift"] = "(theta + 1/2)/x - 2*alpha*x"
    assert verify_construction(changed)["input_sha256"] != report["input_sha256"]


def test_existing_cli_exposes_optional_calculation(tmp_path):
    from fieldbridge.cli import main
    assert main(["verify-construction", str(ROOT / "examples/construction/ito_square.json"),
                 "--out-dir", str(tmp_path)]) == 0


def test_module_command_runs_from_repository(tmp_path):
    result = subprocess.run([sys.executable, "-B", "-m", "fieldbridge", "verify-construction",
                             "examples/construction/ito_square.json", "--out-dir", str(tmp_path)],
                            cwd=ROOT, text=True, capture_output=True, check=True)
    assert json.loads(result.stdout)["candidate"]["passes_local_identity"] is False
