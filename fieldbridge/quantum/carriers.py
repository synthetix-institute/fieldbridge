"""Carriers of quantum realizations and the operators written on them.

A carrier is a Hilbert space together with named operators (tokens). An operator is written as a sum of products of
tokens, for example "X0 X1 + Y0 Y1" or "ad b"; products are applied in the written order (A B means the matrix
product A @ B).

  qubits    {"kind": "qubits", "n": 1..6}                 tokens X0, Y0, Z0, X1, ... and I
  spin      {"kind": "spin", "j": 1/2, 1, 3/2, ... 7/2}   tokens Jx, Jy, Jz, Jp (J+), Jm (J-) and I
  bosons    {"kind": "bosons", "modes": ["a", "b"], "max_quanta": N}
            tokens a, ad (creation), na (number) for every mode name, and I. Each mode is truncated at N quanta;
            within the sector of total number N (declare it as the realization's sector) the truncation is exact for
            products that conserve the number.
  fermions  {"kind": "fermions", "modes": 1..4}           tokens c0, cd0, n0, c1, ... by the Jordan-Wigner
            construction, and I
"""
from __future__ import annotations

from fractions import Fraction
from typing import Dict, List

import numpy as np

PAULI = {"I": np.eye(2, dtype=complex), "X": np.array([[0, 1], [1, 0]], dtype=complex),
         "Y": np.array([[0, -1j], [1j, 0]], dtype=complex), "Z": np.diag([1.0, -1.0]).astype(complex)}


class CarrierError(ValueError):
    pass


def kron_all(mats: List[np.ndarray]) -> np.ndarray:
    out = np.eye(1, dtype=complex)
    for m in mats:
        out = np.kron(out, m)
    return out


def spin_matrices(j: float) -> Dict[str, np.ndarray]:
    """Jx, Jy, Jz, J+ and J- for spin j in the basis m = j, j - 1, ..., -j."""
    m = np.arange(j, -j - 1, -1)
    d = m.size
    jp = np.zeros((d, d), dtype=complex)
    for k in range(1, d):  # J+ |m> = sqrt(j(j+1) - m(m+1)) |m+1>
        jp[k - 1, k] = np.sqrt(j * (j + 1) - m[k] * (m[k] + 1))
    jm = jp.conj().T
    return {"Jx": (jp + jm) / 2, "Jy": (jp - jm) / 2j, "Jz": np.diag(m).astype(complex), "Jp": jp, "Jm": jm}


class Carrier:
    def __init__(self, kind: str, dim: int, tokens: Dict[str, np.ndarray], description: str, spec: Dict):
        self.kind, self.dim, self.tokens, self.description, self.spec = kind, dim, tokens, description, spec

    def product(self, text: str) -> np.ndarray:
        out = np.eye(self.dim, dtype=complex)
        for tok in text.split():
            if tok not in self.tokens:
                raise CarrierError(f"unknown operator '{tok}' on a {self.kind} carrier; known: "
                                   + ", ".join(sorted(self.tokens)))
            out = out @ self.tokens[tok]
        return out

    def operator(self, text: str) -> np.ndarray:
        """A sum of products of tokens, e.g. 'X0 X1 + Y0 Y1'."""
        parts = [p.strip() for p in str(text).split("+")]
        if not all(parts):
            raise CarrierError(f"empty product in operator '{text}'")
        return sum((self.product(p) for p in parts), np.zeros((self.dim, self.dim), dtype=complex))


def make(spec: Dict) -> Carrier:
    kind = spec.get("kind")
    if kind == "qubits":
        n = int(spec.get("n", 0))
        if not 1 <= n <= 6:
            raise CarrierError("qubits: n must be between 1 and 6")
        tokens = {"I": np.eye(2 ** n, dtype=complex)}
        for k in range(n):
            for a in "XYZ":
                tokens[f"{a}{k}"] = kron_all([PAULI[a] if i == k else PAULI["I"] for i in range(n)])
        return Carrier(kind, 2 ** n, tokens, f"{n} spin{'s' if n > 1 else ''}-1/2", spec)
    if kind == "spin":
        j = float(Fraction(str(spec.get("j", "0"))))
        if j <= 0 or j > 3.5 or abs(2 * j - round(2 * j)) > 1e-12:
            raise CarrierError("spin: j must be a positive multiple of 1/2, at most 7/2")
        tokens = spin_matrices(j)
        tokens["I"] = np.eye(int(round(2 * j + 1)), dtype=complex)
        return Carrier(kind, int(round(2 * j + 1)), tokens, f"one spin j = {Fraction(j).limit_denominator()}", spec)
    if kind == "bosons":
        modes = list(spec.get("modes", []))
        nmax = int(spec.get("max_quanta", 0))
        if not modes or len(modes) > 3 or not 1 <= nmax <= 8 or (nmax + 1) ** len(modes) > 729:
            raise CarrierError("bosons: 1-3 modes and 1-8 quanta per mode")
        a1 = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1).astype(complex)
        eye = np.eye(nmax + 1, dtype=complex)
        tokens = {"I": np.eye((nmax + 1) ** len(modes), dtype=complex)}
        for k, name in enumerate(modes):
            a = kron_all([a1 if i == k else eye for i in range(len(modes))])
            tokens[name], tokens[name + "d"], tokens["n" + name] = a, a.conj().T, a.conj().T @ a
        return Carrier(kind, (nmax + 1) ** len(modes), tokens,
                       f"{len(modes)} bosonic modes, at most {nmax} quanta each", spec)
    if kind == "fermions":
        n = int(spec.get("modes", 0))
        if not 1 <= n <= 4:
            raise CarrierError("fermions: 1 to 4 modes")
        lower = np.array([[0, 1], [0, 0]], dtype=complex)  # c|1> = |0> in the basis (|0>, |1>) = (empty, occupied)
        tokens = {"I": np.eye(2 ** n, dtype=complex)}
        for k in range(n):
            c = kron_all([PAULI["Z"] if i < k else lower if i == k else PAULI["I"] for i in range(n)])
            tokens[f"c{k}"], tokens[f"cd{k}"], tokens[f"n{k}"] = c, c.conj().T, c.conj().T @ c
        return Carrier(kind, 2 ** n, tokens, f"{n} fermionic modes", spec)
    raise CarrierError(f"unknown carrier kind '{kind}' (qubits, spin, bosons, fermions)")


def pauli_description(M: np.ndarray, n: int, tol: float = 1e-9) -> List[List]:
    """A qubit operator as a sum of Pauli strings: [[coefficient, 'X0 Z1'], ...]."""
    out = []
    labels = ["I", "X", "Y", "Z"]

    def strings(k):
        if k == 0:
            yield []
            return
        for rest in strings(k - 1):
            for a in labels:
                yield rest + [a]

    for s in strings(n):
        P = kron_all([PAULI[a] for a in s])
        c = np.trace(P @ M) / 2 ** n
        if abs(c) > tol:
            name = " ".join(f"{a}{k}" for k, a in enumerate(s) if a != "I") or "I"
            out.append([round(float(c.real), 12) if abs(c.imag) < tol else complex(c), name])
    return out


def pauli_span(basis: List[np.ndarray], n: int, tol: float = 1e-9) -> List[str]:
    """The Pauli strings that span the same space as the basis, when there are exactly as many as basis elements;
    otherwise an empty list."""
    names: List[str] = []
    for B in basis:
        for c, name in pauli_description(B, n, tol):
            if name not in names and name != "I":
                names.append(name)
    return sorted(names, key=lambda x: (len(x), x)) if len(names) == len(basis) else []
