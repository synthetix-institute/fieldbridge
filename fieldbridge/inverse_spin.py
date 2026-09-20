"""Pair-interaction design with exact constraints and independent dynamics.

Commuting with each variable exchange bond is a sufficient construction route,
not a necessary condition for every possible polarization zero. The response
uses the known collective-phase construction, not an inferred physical law.
"""
from __future__ import annotations

import hashlib
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import sympy as sp


def pair_matrix(sites, i, j, axis):
    pauli = {"I": sp.eye(2), "X": sp.Matrix([[0, 1], [1, 0]]),
             "Y": sp.Matrix([[0, -sp.I], [sp.I, 0]]), "Z": sp.diag(1, -1)}
    return sp.kronecker_product(*(pauli[axis if k in (i, j) else "I"] for k in range(sites)))


def solve_couplings(sites, zz_range="all_pairs"):
    if type(sites) is not int or not 2 <= sites <= 5:
        raise ValueError("sites must be an integer from 2 to 5")
    if zz_range not in ("all_pairs", "nearest_neighbor"):
        raise ValueError("zz_range must be all_pairs or nearest_neighbor")
    pairs = [p for p in combinations(range(sites), 2)
             if zz_range == "all_pairs" or p[1] == p[0]+1]
    terms = [pair_matrix(sites, *p, "Z") for p in pairs]
    exchanges = [pair_matrix(sites, j, j+1, "X")+pair_matrix(sites, j, j+1, "Y")
                 for j in range(sites-1)]
    # Independent exchange variations require a separate commutator per bond.
    blocks = [sp.Matrix.hstack(*[(E*Q-Q*E).reshape(4**sites, 1) for Q in terms])
              for E in exchanges]
    constraints = sp.Matrix.vstack(*blocks)
    basis = constraints.nullspace()
    verified = all(constraints*v == sp.zeros(constraints.rows, 1) for v in basis)
    return {"sites": sites, "zz_range": zz_range, "pairs": [list(p) for p in pairs],
            "coefficient_names": [f"q_{i}{j}" for i, j in pairs],
            "constraint_rank": constraints.rank(), "solution_dimension": len(basis),
            "basis": [[str(x) for x in v] for v in basis], "identities_verified": verified,
            "status": "nonzero_commuting_family" if basis else "no_nonzero_commuting_pair_interaction",
            "requirement": "[X_j X_(j+1)+Y_j Y_(j+1), Q] = 0 for each exchange bond",
            "scope": "Pairwise ZZ ansatz, open XY chain, independently variable exchange strengths."}


def _number(value, name):
    if type(value) not in (int, float, str):
        raise ValueError(f"{name} must be a finite number or rational string")
    try:
        result = sp.Rational(str(value))
    except (ValueError, TypeError, ZeroDivisionError) as error:
        raise ValueError(f"{name} must be finite and rational") from error
    if not result.is_finite:
        raise ValueError(f"{name} must be finite")
    return result


def full_hamiltonian_signal(bonds, zz_coefficients, times):
    """Direct trace for rho=(I+X_0)/2**N, using the full Hamiltonian; hbar=1."""
    sites = len(bonds)+1
    H = sp.zeros(2**sites)
    for j, J in enumerate(bonds):
        H += J*(pair_matrix(sites, j, j+1, "X")+pair_matrix(sites, j, j+1, "Y"))
    for (i, j), q in zz_coefficients.items():
        H += q*pair_matrix(sites, i, j, "Z")
    O = np.kron(np.array([[0, 1], [1, 0]]), np.eye(2**(sites-1)))
    energies, vectors = np.linalg.eigh(np.array(H, dtype=complex))
    weights = np.abs(vectors.conj().T @ O @ vectors)**2 / 2**sites
    gaps = energies[:, None]-energies[None, :]
    return np.array([np.sum(weights*np.cos(gaps*t)).real for t in times])


def collective_prediction(bonds, coupling, times):
    sites = len(bonds)+1
    h = np.zeros((sites, sites))
    for j, J in enumerate(bonds):
        h[j, j+1] = h[j+1, j] = 2*float(J)
    energies, modes = np.linalg.eigh(h)
    returned = np.cos(np.outer(times, energies)) @ np.abs(modes[0])**2
    return returned*np.cos(2*float(coupling)*times)**(sites-1)


def design(spec):
    if not isinstance(spec, dict) or spec.get("schema") != "fieldbridge-spin-design/1":
        raise ValueError("Expected fieldbridge-spin-design/1")
    allowed = {"schema", "question", "bonds", "collective_coupling", "coupling_error", "provenance"}
    if set(spec) != allowed or not spec["question"] or not isinstance(spec["provenance"], dict):
        raise ValueError("Supply question, bonds, collective_coupling, coupling_error and provenance")
    if not isinstance(spec["bonds"], list) or not 2 <= len(spec["bonds"]) <= 4:
        raise ValueError("Supply two to four exchange bonds (three to five spins)")
    bonds = [_number(x, "bond") for x in spec["bonds"]]
    coupling = _number(spec["collective_coupling"], "collective_coupling")
    error = _number(spec["coupling_error"], "coupling_error")
    if coupling == 0 or error == 0:
        raise ValueError("collective_coupling and coupling_error must be nonzero")
    sites = len(bonds)+1
    solution = solve_couplings(sites)
    local = solve_couplings(sites, "nearest_neighbor")
    if solution["basis"] != [["1"]*len(solution["pairs"])]:
        raise ValueError("Solved interaction is not the uniform collective family")
    zero = np.pi/(4*abs(float(coupling)))
    times = np.linspace(0, 2*zero, 129)
    q = {pair: coupling for pair in combinations(range(sites), 2)}
    exact = full_hamiltonian_signal(bonds, q, times)
    predicted = collective_prediction(bonds, coupling, times)
    perturbed = dict(q)
    perturbed[(0, 1)] += error
    broken = float(full_hamiltonian_signal(bonds, perturbed, [zero])[0])
    changed = [J*sp.Rational(j+2, j+1) for j, J in enumerate(bonds)]
    zeros = [float(full_hamiltonian_signal(b, q, [zero])[0]) for b in (bonds, changed)]
    discrepancy = float(np.max(np.abs(exact-predicted)))
    if discrepancy > 1e-10 or max(map(abs, zeros)) > 1e-10:
        raise ValueError("Independent dynamics did not reproduce the prediction")
    return {"schema": "fieldbridge-spin-design-result/1", "status": "calculated",
            "question": spec["question"], "input_sha256": hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest(),
            "provenance": spec["provenance"], "novelty_established": False, "source_alignment_verified": False,
            "units": "hbar=1; all couplings in the same energy unit", "construction": solution,
            "nearest_neighbor_comparison": local,
            "prediction": {"preparation": "rho(0)=(I+X_0)/2**N", "observable": "X_0",
                "expression": "Re[exp(-i h t)]_00 * cos(2 lambda t)**(N-1); h_(j,j+1)=2 J_j",
                "first_collective_zero_time": float(zero), "first_collective_zero_formula": "pi/(4*abs(lambda))",
                "scope": "End-spin polarization, stated preparation, uniform all-pair ZZ coupling; earlier exchange zeros may occur."},
            "checks": {"max_full_hamiltonian_error": discrepancy, "exchange_disorder_zero_values": zeros,
                "changed_exchange_bonds": [str(J) for J in changed], "single_pair_coupling_error": str(error),
                "perturbed_signal_at_original_zero": broken, "uniform_coupling_control_fires": bool(abs(broken)>1e-8)},
            "response": [{"time": float(t), "predicted": float(p), "full_hamiltonian": float(v)}
                         for t, p, v in zip(times, predicted, exact)]}


def render(report):
    local, check = report["nearest_neighbor_comparison"], report["checks"]
    return ("# Designing a spin cancellation\n\n"+report["question"]+"\n\n"
        "The exact constraints give one free coefficient: Q=lambda sum_(i<j) Z_i Z_j.\n\n"
        f"Nearest-neighbour Ising interactions leave {local['solution_dimension']} free coefficients. "
        "Zero means no nonzero solution within this commuting pair-interaction construction, not impossibility of every cancellation mechanism.\n\n"
        "For rho(0)=(I+X_0)/2**N the polarization is Re[exp(-i h t)]_00 cos(2 lambda t)**(N-1), with h_(j,j+1)=2 J_j and hbar=1.\n\n"
        f"First collective zero: t={report['prediction']['first_collective_zero_time']:.8g}. "
        f"Full-Hamiltonian error across 129 times: {check['max_full_hamiltonian_error']:.3g}.\n\n"
        f"Varying exchange bonds gives signals {check['exchange_disorder_zero_values']} at that time. "
        f"Changing one Ising coupling by {check['single_pair_coupling_error']} gives {check['perturbed_signal_at_original_zero']:.8g}. "
        f"Uniform-coupling control fires: {check['uniform_coupling_control_fires']}.\n\n"
        "The coupling constraint is solved. The response uses an established analytic construction checked by direct Hamiltonian evolution. "
        "This authored example makes no claim of literature extraction, experiment or established novelty.\n")


def run_file(path, out_dir):
    destination = Path(out_dir)
    destination.mkdir(parents=True, exist_ok=True)
    report_path = destination / "design.json"
    report_path.write_text(json.dumps({"status": "running"})+"\n", encoding="utf-8")
    (destination / "design.md").write_text("Calculation running.\n", encoding="utf-8")
    try:
        spec = json.loads(Path(path).read_text(encoding="utf-8"))
        report = design(spec)
        report["software"] = {"numpy": np.__version__, "sympy": sp.__version__}
        report["implementation_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        (destination / "input.json").write_text(json.dumps(spec, indent=2)+"\n", encoding="utf-8")
        report_path.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
        (destination / "design.md").write_text(render(report), encoding="utf-8")
        return report
    except Exception as error:
        report_path.write_text(json.dumps({"status": "failed", "reason": str(error)})+"\n", encoding="utf-8")
        (destination / "design.md").write_text("Calculation failed: "+str(error)+"\n", encoding="utf-8")
        (destination / "input.json").write_text("null\n", encoding="utf-8")
        raise
