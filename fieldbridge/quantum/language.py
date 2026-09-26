"""The language of mechanisms on quantum carriers: identity, detachment, attachment and co-discovery.

A realization is I = ((Omega, Xi); C, R, P; A):

  Xi     the carrier: a Hilbert space and the operators written on it (carriers.py), or one symmetry sector of it
  Omega  the operation: the Hamiltonian, a sum of terms
  C      the closure, what is specified or discarded to close the equations: closed unitary dynamics, hbar = 1;
         a declared sector keeps a conserved quantity fixed
  R      the observable
  P      the protocol: the state is prepared in the top eigenstate of the observable and H acts from t = 0
  A      the parameters, through which a material or apparatus implements the model

Target: the Bloch rotation. The Hamiltonian and the observable generate the Lie algebra su(2),
[J_a, J_b] = i eps_abc J_c. The expectation values m_a = <J_a> then rotate, dm/dt = Omega x m, and the observable
follows the Rabi law f(t) = cos^2 theta + sin^2 theta cos(|Omega| t), with theta the angle between Omega and the
observable, on every carrier.

Letters of a derivation:

  S  sector      a conserved operator restricts the carrier to one of its eigenspaces
  A  algebra     commutators of the Hamiltonian and the observable are added until their span closes: the dynamical
                 Lie algebra
  K  canonical   a basis J_1, J_2, J_3 with [J_a, J_b] = i eps_abc J_c is found from the Killing form, with J_3 along
                 Omega and J_1 in the plane of Omega and the observable
  L  law         the observable is evolved exactly on the carrier and compared with the Rabi law

Detachment keeps the part of a realization that does not depend on the carrier: the algebra, |Omega|, the weight
|n| of the observable and the angle theta. Attachment writes that mechanism on another carrier in the carrier's own
operators (a tunnelling, an exchange coupling, a pairing field) and derives it there again.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from math import acos, cos, degrees, pi, sin, sqrt
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from .carriers import Carrier, CarrierError, make, pauli_description, pauli_span

SCHEMA = "fieldbridge-quantum/1"
REQUIRED = ("name", "question", "assumptions", "provenance", "carrier", "hamiltonian", "observable")
LETTERS = {"S": "sector", "A": "algebra", "K": "canonical form", "L": "law"}
SLOTS = {"S": ("Xi, C", "a conserved operator restricts the carrier to one of its eigenspaces"),
         "A": ("Omega, R", "commutators of the Hamiltonian and the observable are added until the span closes"),
         "K": ("Omega, R", "a basis with [J_a, J_b] = i eps_abc J_c is found; J_3 along Omega"),
         "L": ("P, R", "the observable is evolved exactly on the carrier and compared with the Rabi law")}
MAX_DIM = 40
TOL = 1e-9


class SpecError(ValueError):
    pass


@dataclass
class Realization:
    name: str
    field: str
    spec: Dict
    carrier: Carrier
    params: Dict[str, float]
    terms: List[Tuple[str, np.ndarray]]
    H: np.ndarray
    O: np.ndarray
    sector: Optional[Dict]


# ------------------------------------------------------------------------------------------------ specifications
def _coefficient(value, params: Dict[str, float]) -> float:
    if isinstance(value, (int, float)):
        return float(value)
    import sympy as sp
    symbols = {k: sp.Symbol(k) for k in params}
    try:
        expr = sp.sympify(str(value), locals=symbols)
    except (sp.SympifyError, TypeError) as err:
        raise SpecError(f"cannot read the coefficient '{value}'") from err
    unknown = {str(s) for s in expr.free_symbols} - set(params)
    if unknown:
        raise SpecError(f"coefficient '{value}' uses undeclared parameters: {', '.join(sorted(unknown))}")
    val = complex(expr.subs({symbols[k]: v for k, v in params.items()}).evalf())
    if abs(val.imag) > 1e-12:
        raise SpecError(f"coefficient '{value}' is not real")
    return float(val.real)


def _hermitian(M: np.ndarray) -> bool:
    return bool(np.allclose(M, M.conj().T, atol=1e-10 * max(1.0, float(np.abs(M).max()))))


def _terms(entries, carrier: Carrier, params, what: str) -> List[Tuple[str, np.ndarray]]:
    if not isinstance(entries, list) or not entries:
        raise SpecError(f"{what} must be a non-empty list of terms")
    out = []
    for e in entries:
        if not isinstance(e, dict) or "operator" not in e or "coefficient" not in e:
            raise SpecError(f"every {what} term needs 'coefficient' and 'operator'")
        c = _coefficient(e["coefficient"], params)
        try:
            M = carrier.operator(e["operator"])
        except CarrierError as err:
            raise SpecError(str(err)) from err
        if e.get("hc"):
            M = M + M.conj().T
        out.append((f"{e['coefficient']} [{e['operator']}{' + h.c.' if e.get('hc') else ''}]", c * M))
    if not _hermitian(sum(m for _, m in out)):
        raise SpecError(f"the {what} is not Hermitian; mark non-Hermitian products with \"hc\": true")
    return out


def load(source: Union[str, Path, Dict]) -> Realization:
    spec = source if isinstance(source, dict) else json.loads(Path(source).read_text(encoding="utf-8"))
    if spec.get("schema") != SCHEMA:
        raise SpecError(f"schema must be {SCHEMA}")
    missing = [k for k in REQUIRED if not spec.get(k)]
    if missing:
        raise SpecError("missing: " + ", ".join(missing))
    try:
        carrier = make(spec["carrier"])
    except CarrierError as err:
        raise SpecError(str(err)) from err
    params = {k: float(v) for k, v in (spec.get("parameters") or {}).items()}
    terms = _terms(spec["hamiltonian"], carrier, params, "Hamiltonian")
    O = sum(m for _, m in _terms(spec["observable"], carrier, params, "observable"))
    sector = None
    if spec.get("sector"):
        sec = spec["sector"]
        S = sum(m for _, m in _terms(sec["operator"], carrier, params, "sector operator"))
        sector = {"operator": " + ".join(f"{e['coefficient']} [{e['operator']}]" for e in sec["operator"]),
                  "value": _coefficient(sec["value"], params), "matrix": S}
    return Realization(name=spec["name"], field=str(spec.get("field", "unspecified")), spec=spec, carrier=carrier,
                       params=params, terms=terms, H=sum(m for _, m in terms), O=O, sector=sector)


# ------------------------------------------------------------------------------------------------ the letters
def restrict(real: Realization, H: Optional[np.ndarray] = None):
    """S: the eigenspace of the sector operator with the declared value. Returns the isometry V (columns span the
    sector) and a record; V is the identity when no sector is declared."""
    H = real.H if H is None else H
    d = real.carrier.dim
    if real.sector is None:
        return np.eye(d, dtype=complex), None
    S = real.sector["matrix"]
    scale = max(1.0, float(np.abs(H).max()), float(np.abs(real.O).max()))
    for name, M in (("Hamiltonian", H), ("observable", real.O)):
        if float(np.abs(M @ S - S @ M).max()) > 1e-9 * scale * max(1.0, float(np.abs(S).max())):
            raise SpecError(f"the sector operator does not commute with the {name}: it is not conserved")
    w, U = np.linalg.eigh(S)
    V = U[:, np.abs(w - real.sector["value"]) < 1e-8]
    if V.shape[1] == 0:
        raise SpecError(f"{real.sector['value']} is not an eigenvalue of the sector operator")
    return V, {"operator": real.sector["operator"], "value": real.sector["value"], "dimension": int(V.shape[1]),
               "full_dimension": d}


def _inner(A: np.ndarray, B: np.ndarray) -> float:
    return float(np.real(np.trace(A @ B)))


def _traceless(A: np.ndarray) -> np.ndarray:
    return A - np.trace(A) / A.shape[0] * np.eye(A.shape[0])


def bracket(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """-i[A, B]: Hermitian for Hermitian A and B."""
    return -1j * (A @ B - B @ A)


def lie_closure(gens: List[np.ndarray], max_dim: int = MAX_DIM) -> Tuple[List[np.ndarray], bool]:
    """A: an orthonormal basis (trace inner product) of the real Lie algebra generated by the traceless parts of the
    generators under -i[ , ]. Returns the basis and whether it closed within max_dim."""
    basis: List[np.ndarray] = []

    def add(M: np.ndarray) -> None:
        v = M.copy()
        for _ in range(2):  # Gram-Schmidt twice for stability
            for B in basis:
                v = v - _inner(B, v) * B
        nv = sqrt(max(_inner(v, v), 0.0))
        if nv > TOL * max(1.0, sqrt(max(_inner(M, M), 0.0))):
            basis.append((v + v.conj().T) / (2 * nv))

    for g in gens:
        add(_traceless(g))
    i = 0
    while i < len(basis):
        for j in range(i):
            add(bracket(basis[i], basis[j]))
            if len(basis) > max_dim:
                return basis, False
        i += 1
    return basis, True


def structure_constants(basis: List[np.ndarray]) -> np.ndarray:
    d = len(basis)
    f = np.zeros((d, d, d))
    for a in range(d):
        for b in range(d):
            C = bracket(basis[a], basis[b])
            f[a, b] = [_inner(C, basis[c]) for c in range(d)]
    return f


def killing(f: np.ndarray) -> np.ndarray:
    """K_ab = tr(ad_a ad_b), with (ad_a)_cd = f_adc."""
    ad = np.transpose(f, (0, 2, 1))  # ad[a][c, d] = f[a, d, c]
    return np.einsum("acd,bdc->ab", ad, ad)


def canonical_su2(basis: List[np.ndarray], H0: np.ndarray, O0: np.ndarray) -> Optional[Dict[str, object]]:
    """K: for a three-dimensional algebra with negative definite Killing form, the basis J with
    [J_a, J_b] = i eps_abc J_c in the frame J_3 along Omega, J_1 in the plane of Omega and the observable."""
    if len(basis) != 3:
        return None
    f = structure_constants(basis)
    K = killing(f)
    w, U = np.linalg.eigh(K)
    if not np.all(w < -1e-9 * max(1.0, float(np.abs(w).max()))):
        return None
    M = U @ np.diag(np.sqrt(-w / 2)) @ U.T  # e = M J with M M^T = -K/2
    Minv = np.linalg.inv(M)
    J = [sum(Minv[a, b] * basis[b] for b in range(3)) for a in range(3)]
    if _inner(bracket(J[0], J[1]), J[2]) < 0:  # orientation: make the structure constants +eps
        J = [-X for X in J]
    norm = [_inner(X, X) for X in J]
    Om = np.array([_inner(H0, J[a]) / norm[a] for a in range(3)])
    nv = np.array([_inner(O0, J[a]) / norm[a] for a in range(3)])
    W, N = float(np.linalg.norm(Om)), float(np.linalg.norm(nv))
    if W < 1e-12 or N < 1e-12:
        return None
    e3 = Om / W
    perp = nv - (nv @ e3) * e3
    if np.linalg.norm(perp) < 1e-12 * N:
        perp = np.cross(e3, [1.0, 0.0, 0.0]) if abs(e3[0]) < 0.9 else np.cross(e3, [0.0, 1.0, 0.0])
    e1 = perp / np.linalg.norm(perp)
    e2 = np.cross(e3, e1)
    R = np.array([e1, e2, e3])
    Jc = [sum(R[a, b] * J[b] for b in range(3)) for a in range(3)]
    theta = acos(max(-1.0, min(1.0, float(nv @ Om) / (N * W))))
    residual_H = float(np.linalg.norm(H0 - sum(Om[a] * J[a] for a in range(3))))
    residual_O = float(np.linalg.norm(O0 - sum(nv[a] * J[a] for a in range(3))))
    cas = sum(X @ X for X in Jc)
    lam = np.linalg.eigvalsh((cas + cas.conj().T) / 2)
    js = np.round(2 * (-1 + np.sqrt(1 + 4 * np.maximum(lam, 0))) / 2) / 2
    content = {}
    for jv in js:
        content[float(jv)] = content.get(float(jv), 0) + 1
    rep = [{"j": j, "copies": int(round(c / (2 * j + 1)))} for j, c in sorted(content.items(), reverse=True)]
    return {"J": Jc, "rate": W, "weight": N, "theta": theta, "structure_residual": _structure_residual(Jc),
            "decomposition_residual": max(residual_H, residual_O), "representation": rep}


def _structure_residual(J: List[np.ndarray]) -> float:
    eps = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1}
    worst = 0.0
    for (a, b, c) in eps:
        worst = max(worst, float(np.abs(J[a] @ J[b] - J[b] @ J[a] - 1j * J[c]).max()))
    return worst


def rabi_law(H: np.ndarray, O: np.ndarray, canon: Dict, n_times: int = 161) -> Dict[str, object]:
    """L: exact evolution from the top eigenstate of the observable, against
    <O>(t) = o0 + |n| j (cos^2 theta + sin^2 theta cos(|Omega| t))."""
    o0 = float(np.real(np.trace(O))) / O.shape[0]
    w, U = np.linalg.eigh(O)
    psi0 = U[:, -1]
    j_top = (float(w[-1]) - o0) / canon["weight"]
    E, V = np.linalg.eigh(H)
    c = V.conj().T @ psi0
    W, th = canon["rate"], canon["theta"]
    t = np.linspace(0.0, 4 * pi / W, n_times)
    exact = np.empty_like(t)
    for k, tk in enumerate(t):
        psi = V @ (np.exp(-1j * E * tk) * c)
        exact[k] = float(np.real(psi.conj() @ O @ psi))
    f_exact = (exact - o0) / (canon["weight"] * j_top)
    f_law = cos(th) ** 2 + sin(th) ** 2 * np.cos(W * t)
    return {"times": t.tolist(), "f_exact": f_exact.tolist(), "f_law": f_law.tolist(), "j_top": j_top,
            "residual": float(np.max(np.abs(f_exact - f_law))), "inversion_time": pi / W,
            "value_at_inversion": cos(th) ** 2 - sin(th) ** 2}


def _describe_algebra(dim: int, sector_dim: int) -> str:
    if dim == 1:
        return "the observable commutes with the Hamiltonian: it is conserved and nothing rotates"
    if dim == sector_dim ** 2 - 1 and sector_dim > 2:
        return (f"the algebra is su({sector_dim}), dimension {dim}: every observable of the {sector_dim}-dimensional "
                f"carrier is reached, not only a rotating three-vector")
    return f"the algebra has dimension {dim}, not 3: the observable is not carried by one rotating three-vector"


def derive_bloch_rotation(real: Realization, law: bool = True) -> Dict[str, object]:
    """Derive the Bloch rotation in one realization; returns the word, the steps, the signature and the law."""
    word: List[str] = []
    steps: List[Dict[str, object]] = []
    out: Dict[str, object] = {"name": real.name, "field": real.field, "carrier": real.carrier.description,
                              "status": "obstructed"}
    V, sec = restrict(real)
    if sec is not None:
        word.append("S")
        steps.append({"letter": "S", **sec})
    Hs, Os = V.conj().T @ real.H @ V, V.conj().T @ real.O @ V
    basis, closed = lie_closure([Hs, Os])
    word.append("A")
    step_a = {"letter": "A", "dimension": len(basis) if closed else f"> {MAX_DIM}", "closed": closed}
    if real.carrier.kind == "qubits" and sec is None and closed and len(basis) <= 15:
        step_a["closing_operators"] = pauli_span(basis, real.carrier.spec["n"])
    steps.append(step_a)
    canon = canonical_su2(basis, _traceless(Hs), _traceless(Os)) if closed else None
    if canon is None:
        dim = len(basis) if closed else MAX_DIM + 1
        reason = _describe_algebra(dim, Hs.shape[0]) if closed else f"the algebra exceeds dimension {MAX_DIM}"
        if dim == 3:
            reason = "the three-dimensional algebra is not su(2) (its Killing form is not negative definite)"
        cause = _single_term_cause(real, V) if dim > 3 and len(real.terms) > 1 else None
        if cause:
            reason += f"; without the term {cause} the algebra is su(2): that term is the obstruction"
        out.update(obstruction=reason, algebra_dimension=dim)
        return {**out, "word": "".join(word), "class": " ".join(word), "steps": steps}
    word.append("K")
    steps.append({"letter": "K", "rate": canon["rate"], "weight": canon["weight"], "theta_deg": degrees(canon["theta"]),
                  "representation": canon["representation"], "structure_residual": canon["structure_residual"],
                  "decomposition_residual": canon["decomposition_residual"],
                  **({"generators": [pauli_description(V @ X @ V.conj().T, real.carrier.spec["n"]) for X in canon["J"]]}
                     if real.carrier.kind == "qubits" and sec is None else {})})
    out.update(status="reached", algebra_dimension=3)
    lawr = {}
    if law:
        lawr = rabi_law(Hs, Os, canon)
        word.append("L")
        steps.append({"letter": "L", "residual": lawr["residual"], "j_top": lawr["j_top"],
                      "inversion_time": lawr["inversion_time"], "value_at_inversion": lawr["value_at_inversion"]})
    signature = {"algebra": "su(2)", "relations": "[J_a, J_b] = i eps_abc J_c", "rate": canon["rate"],
                 "weight": canon["weight"], "theta_deg": degrees(canon["theta"]),
                 "law": "f(t) = cos^2 theta + sin^2 theta cos(|Omega| t)"}
    return {**out, "word": "".join(word), "class": " ".join(word), "steps": steps, "signature": signature,
            "representation": canon["representation"], "sector_dimension": int(Hs.shape[0]), "law": lawr}


def _single_term_cause(real: Realization, V: np.ndarray) -> Optional[str]:
    """The Hamiltonian term whose removal leaves su(2), if there is exactly such a term."""
    found = []
    Os = V.conj().T @ real.O @ V
    for k, (label, M) in enumerate(real.terms):
        rest = sum((m for i, (_, m) in enumerate(real.terms) if i != k), np.zeros_like(real.H))
        Hs = V.conj().T @ rest @ V
        basis, closed = lie_closure([Hs, Os])
        if closed and canonical_su2(basis, _traceless(Hs), _traceless(Os)) is not None:
            found.append(label)
    return found[0] if len(found) == 1 else None


# ------------------------------------------------------------------------------------------------ detach / attach
def detach(real: Realization) -> Dict[str, object]:
    """The carrier-free part of a realization (the signature) and the carrier-specific part that is left behind."""
    row = derive_bloch_rotation(real, law=True)
    if row["status"] != "reached":
        return {"name": real.name, "detached": False, "reason": row["obstruction"], "derivation": row}
    k = next(s for s in row["steps"] if s["letter"] == "K")
    a = next(s for s in row["steps"] if s["letter"] == "A")
    return {"name": real.name, "detached": True, "signature": row["signature"],
            "left_behind": {"carrier": real.carrier.description, "field": real.field,
                            "hilbert_dimension": real.carrier.dim, "sector_dimension": row["sector_dimension"],
                            "representation": row["representation"],
                            **({"closing_operators": a["closing_operators"]} if a.get("closing_operators") else {}),
                            **({"generators": k["generators"]} if "generators" in k else {})},
            "law": {"residual": row["law"]["residual"], "inversion_time": row["law"]["inversion_time"]},
            "derivation": row}


CARRIERS = {
    "qubit": "one spin-1/2 in a field (nuclear magnetic resonance)",
    "spin": "one spin j (size = 2j)",
    "bosons": "two bosonic modes with N particles (size = N; a Bose-Josephson junction)",
    "chain": "an exchange chain of N spins with one flipped spin (size = N; state transfer)",
    "fermion-pair": "two fermionic modes k and -k (a Cooper-pair level; Anderson pseudospin)",
    "correlated-pair": "two spins whose correlations form the rotating vector (the Module 11 spins)",
    "collective": "N spins-1/2 driven together (size = N; collective spin)",
}
FIELDS = {"qubit": "nuclear magnetic resonance", "spin": "spin physics", "bosons": "cold atoms",
          "chain": "quantum information", "fermion-pair": "superconductivity", "correlated-pair": "quantum information",
          "collective": "spin physics"}


def attached_spec(signature: Dict, to: str, size: Optional[int] = None, source: str = "") -> Dict:
    """Write the detached mechanism on another carrier: H = |Omega| (sin theta D + cos theta Z) and O = |n| Z, where
    Z is the carrier's natural observable and D its natural drive, with [D, T] = i Z for a third generator T."""
    W, N, th = signature["rate"], signature["weight"], signature["theta_deg"] * pi / 180
    Ws, Wc = W * sin(th), W * cos(th)
    if to not in CARRIERS:
        raise SpecError(f"unknown carrier '{to}'; choose from: {', '.join(CARRIERS)}")
    r = lambda x: float(f"{x:.12g}")  # noqa: E731
    base = {"schema": SCHEMA, "field": FIELDS[to],
            "question": f"Does the rotation detached from '{source}' hold on {CARRIERS[to]}?",
            "assumptions": ["Closed unitary dynamics; hbar = 1.",
                            f"The Hamiltonian writes the detached rotation (rate {W:.6g}, angle "
                            f"{signature['theta_deg']:.4g} degrees) in the carrier's own operators."],
            "provenance": {"origin": f"attached from '{source}' by fieldbridge quantum attach",
                           "novelty": "not established; a known representation of su(2)"}}
    if to == "qubit":
        return {**base, "name": f"{source} on one spin-1/2", "carrier": {"kind": "qubits", "n": 1},
                "hamiltonian": [{"coefficient": r(Ws / 2), "operator": "X0"}, {"coefficient": r(Wc / 2), "operator": "Z0"}],
                "observable": [{"coefficient": r(N / 2), "operator": "Z0"}],
                "native": {"rabi_frequency": r(Ws), "detuning": r(Wc)}}
    if to == "spin":
        j = (size or 2) / 2
        return {**base, "name": f"{source} on a spin {j:g}", "carrier": {"kind": "spin", "j": f"{size or 2}/2"},
                "hamiltonian": [{"coefficient": r(Ws), "operator": "Jx"}, {"coefficient": r(Wc), "operator": "Jz"}],
                "observable": [{"coefficient": r(N), "operator": "Jz"}], "native": {"field_x": r(Ws), "field_z": r(Wc)}}
    if to == "bosons":
        n = size or 4
        return {**base, "name": f"{source} on two bosonic modes with {n} particles",
                "carrier": {"kind": "bosons", "modes": ["a", "b"], "max_quanta": n},
                "sector": {"operator": [{"coefficient": 1, "operator": "na + nb"}], "value": n},
                "hamiltonian": [{"coefficient": r(Ws / 2), "operator": "ad b", "hc": True},
                                {"coefficient": r(Wc / 2), "operator": "na"}, {"coefficient": r(-Wc / 2), "operator": "nb"}],
                "observable": [{"coefficient": r(N / 2), "operator": "na"}, {"coefficient": r(-N / 2), "operator": "nb"}],
                "native": {"tunnelling": r(Ws / 2), "energy_difference": r(Wc)}}
    if to == "chain":
        n = size or 4
        S = (n - 1) / 2
        kappa = [sqrt((j + 1) * (n - 1 - j)) / 2 for j in range(n - 1)]
        bonds = [r(Ws * k / 2) for k in kappa]
        ham = [{"coefficient": bonds[j], "operator": f"X{j} X{j + 1} + Y{j} Y{j + 1}"} for j in range(n - 1)]
        ham += [{"coefficient": r(-Wc * (j - S) / 2), "operator": f"Z{j}"} for j in range(n) if abs(j - S) > 1e-12]
        obs = [{"coefficient": r(-N * (j - S) / 2), "operator": f"Z{j}"} for j in range(n) if abs(j - S) > 1e-12]
        return {**base, "name": f"{source} on an exchange chain of {n} spins", "carrier": {"kind": "qubits", "n": n},
                "sector": {"operator": [{"coefficient": 1, "operator": f"Z{j}"} for j in range(n)], "value": n - 2},
                "hamiltonian": ham if Wc else ham[: n - 1], "observable": obs,
                "native": {"exchange_couplings": bonds, "coupling_ratios": [r(k / kappa[0]) for k in kappa],
                           "site_fields": [r(-Wc * (j - S) / 2) for j in range(n)],
                           "observable": "position of the flipped spin, sum_j (j - (N-1)/2) (1 - Z_j)/2"}}
    if to == "fermion-pair":
        return {**base, "name": f"{source} on a Cooper-pair level", "carrier": {"kind": "fermions", "modes": 2},
                "hamiltonian": [{"coefficient": r(Wc / 2), "operator": "n0"}, {"coefficient": r(Wc / 2), "operator": "n1"},
                                {"coefficient": r(Ws / 2), "operator": "cd0 cd1", "hc": True}],
                "observable": [{"coefficient": r(N / 2), "operator": "n0"}, {"coefficient": r(N / 2), "operator": "n1"}],
                "native": {"single_particle_energy": r(Wc / 2), "pairing_amplitude": r(-Ws / 2)}}
    if to == "correlated-pair":
        return {**base, "name": f"{source} on the correlations of two spins", "carrier": {"kind": "qubits", "n": 2},
                "hamiltonian": [{"coefficient": r(Ws / 2), "operator": "Z0 Z1"}, {"coefficient": r(Wc / 2), "operator": "X0"}],
                "observable": [{"coefficient": r(N / 2), "operator": "X0"}],
                "native": {"ising_coupling_g": r(Ws / 2), "transverse_field_h": r(Wc / 2)}}
    n = size or 3
    return {**base, "name": f"{source} on {n} spins driven together", "carrier": {"kind": "qubits", "n": n},
            "hamiltonian": [{"coefficient": r(Ws / 2), "operator": " + ".join(f"X{k}" for k in range(n))},
                            {"coefficient": r(Wc / 2), "operator": " + ".join(f"Z{k}" for k in range(n))}],
            "observable": [{"coefficient": r(N / 2), "operator": " + ".join(f"Z{k}" for k in range(n))}],
            "native": {"field_x": r(Ws / 2), "field_z": r(Wc / 2)}}


def attach(source: Realization, to: str, size: Optional[int] = None) -> Dict[str, object]:
    """Detach the rotation from the source, write it on the carrier 'to', derive it there and compare."""
    d = detach(source)
    if not d["detached"]:
        return {"attached": False, "reason": d["reason"], "source": source.name}
    spec = attached_spec(d["signature"], to, size, source.name)
    target = load(spec)
    row = derive_bloch_rotation(target)
    sig_t = row.get("signature", {})
    preserved = {k: {"source": d["signature"][k], "target": sig_t.get(k)} for k in ("rate", "weight", "theta_deg")}
    return {"attached": row["status"] == "reached", "source": source.name, "to": to, "spec": spec,
            "native": spec["native"], "preserved": preserved,
            "changed": {"carrier": {"source": d["left_behind"]["carrier"], "target": target.carrier.description},
                        "sector_dimension": {"source": d["left_behind"]["sector_dimension"],
                                             "target": row.get("sector_dimension")},
                        "representation": {"source": d["left_behind"]["representation"],
                                           "target": row.get("representation")}},
            "law": {"source_residual": d["law"]["residual"], "target_residual": row.get("law", {}).get("residual"),
                    "inversion_time": row.get("law", {}).get("inversion_time")},
            "source_derivation": d["derivation"], "target_derivation": row}


# ------------------------------------------------------------------------------------------------ co-discovery
def codiscover(reals: List[Realization]) -> Dict[str, object]:
    rows = [derive_bloch_rotation(r) for r in reals]
    reached = [r for r in rows if r["status"] == "reached"]
    classes: Dict[str, List[str]] = {}
    for r in reached:
        classes.setdefault(r["class"], []).append(r["name"])
    summary = {"target": "Bloch rotation (the Hamiltonian and the observable generate su(2))",
               "reached": len(reached), "obstructed": len(rows) - len(reached),
               "fields_reached": sorted({r["field"] for r in reached}), "derivation_classes": classes,
               "law_residual_max": max((r["law"]["residual"] for r in reached), default=None),
               "structure_residual_max": max((next(s for s in r["steps"] if s["letter"] == "K")["structure_residual"]
                                              for r in reached), default=None)}
    return {"summary": summary, "rows": rows}


def _rep(rep: List[Dict]) -> str:
    def j_text(j):
        return f"{int(j)}" if float(j).is_integer() else f"{int(round(2 * j))}/2"
    return " + ".join(f"{c} x j={j_text(r['j'])}" if (c := r["copies"]) > 1 else f"j={j_text(r['j'])}" for r in rep)


def markdown(report: Dict[str, object]) -> str:
    s = report["summary"]
    lines = ["# Co-discovery by construction: the Bloch rotation", "",
             f"Target: {s['target']}. Reached in {s['reached']} realizations from {len(s['fields_reached'])} fields "
             f"({', '.join(s['fields_reached'])}); obstructed in {s['obstructed']}.", "",
             "| realization | field | carrier | derivation | algebra dimension | representation | rate | angle | "
             "law residual |", "|---|---|---|---|---:|---|---:|---:|---:|"]
    for r in report["rows"]:
        if r["status"] == "reached":
            sig = r["signature"]
            lines.append(f"| {r['name']} | {r['field']} | {r['carrier']} | {r['word']} | 3 | {_rep(r['representation'])} | "
                         f"{sig['rate']:.4g} | {sig['theta_deg']:.1f} | {r['law']['residual']:.1e} |")
        else:
            lines.append(f"| {r['name']} | {r['field']} | {r['carrier']} | {r['word']} | {r['algebra_dimension']} | - | "
                         f"- | - | - |")
    lines += ["", "## Obstructions", ""]
    for r in report["rows"]:
        if r["status"] != "reached":
            lines.append(f"- {r['name']} ({r['field']}): {r['obstruction']}")
    lines += ["", "## Invariants", "",
              f"- Algebra: su(2) in every reached realization; largest deviation from [J_a, J_b] = i eps_abc J_c: "
              f"{s['structure_residual_max']:.1e}.",
              f"- Law: the observable follows cos^2 theta + sin^2 theta cos(|Omega| t) on every carrier; largest "
              f"deviation of the exact evolution: {s['law_residual_max']:.1e}.", "",
              "## Letters and the slots they act on", "", "| letter | transformation | slots | action |",
              "|---|---|---|---|"]
    lines += [f"| {k} | {LETTERS[k]} | {SLOTS[k][0]} | {SLOTS[k][1]} |" for k in LETTERS]
    return "\n".join(lines) + "\n"
