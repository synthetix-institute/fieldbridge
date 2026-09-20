import json
from itertools import combinations
from pathlib import Path
import subprocess
import sys

import pytest

np = pytest.importorskip("numpy")
sp = pytest.importorskip("sympy")
from fieldbridge.inverse_spin import collective_prediction, design, full_hamiltonian_signal, solve_couplings

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/construction/spin_cancellation_design.json"


@pytest.mark.parametrize("sites", [2, 3, 4, 5])
def test_exact_inverse_constraints_find_uniform_coupling(sites):
    result = solve_couplings(sites)
    assert result["solution_dimension"] == 1
    assert result["basis"] == [["1"]*(sites*(sites-1)//2)]
    assert result["identities_verified"]
    assert result["constraint_rank"] == sites*(sites-1)//2-1


@pytest.mark.parametrize("sites", [2, 3, 4, 5])
def test_local_restriction_changes_admissibility(sites):
    result = solve_couplings(sites, "nearest_neighbor")
    assert result["solution_dimension"] == (1 if sites == 2 else 0)


@pytest.mark.parametrize("bonds", [[1, "7/10"], [3, -4, 2], [1, 2, 0, 3]])
def test_full_dynamics_independently_reproduces_response(bonds):
    bonds = [sp.Rational(J) for J in bonds]
    coupling = sp.Rational(3, 10)
    times = np.linspace(0, 7, 37)
    q = {pair: coupling for pair in combinations(range(len(bonds)+1), 2)}
    error = np.max(np.abs(full_hamiltonian_signal(bonds, q, times)-collective_prediction(bonds, coupling, times)))
    assert error < 1e-12


def test_worked_example_solves_predicts_and_exposes_a_broken_assumption():
    result = design(json.loads(EXAMPLE.read_text()))
    assert result["checks"]["max_full_hamiltonian_error"] < 1e-12
    assert abs(result["prediction"]["first_collective_zero_time"]-np.pi/5) < 1e-14
    assert result["checks"]["uniform_coupling_control_fires"]
    assert result["checks"]["perturbed_signal_at_original_zero"] == pytest.approx(0.001689, abs=1e-6)
    assert max(map(abs, result["checks"]["exchange_disorder_zero_values"])) < 1e-12
    assert result["nearest_neighbor_comparison"]["solution_dimension"] == 0
    assert result["novelty_established"] is False


@pytest.mark.parametrize("key,value", [("schema", "bad"), ("bonds", [1]), ("bonds", [1, "nan"]),
                                      ("collective_coupling", 0), ("coupling_error", 0), ("unexpected", 1)])
def test_invalid_design_refuses(key, value):
    spec = json.loads(EXAMPLE.read_text())
    spec[key] = value
    with pytest.raises(ValueError):
        design(spec)


def test_cli_rerun_cannot_leave_a_stale_success(tmp_path):
    source = tmp_path / "input.json"
    source.write_text(EXAMPLE.read_text())
    out = tmp_path / "out"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "design-spin-cancellation", str(source), "--out-dir", str(out)]
    run = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    report = json.loads((out / "design.json").read_text())
    assert report["status"] == "calculated"
    assert len(report["implementation_sha256"]) == 64
    assert "First collective zero" in (out / "design.md").read_text()
    source.write_text("{}")
    run = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    assert run.returncode != 0
    assert json.loads((out / "design.json").read_text())["status"] == "failed"
    assert json.loads((out / "input.json").read_text()) is None


def test_solver_resource_and_model_bounds():
    for size in [True, 1, 6, "3"]:
        with pytest.raises(ValueError):
            solve_couplings(size)
    with pytest.raises(ValueError):
        solve_couplings(3, "arbitrary")
