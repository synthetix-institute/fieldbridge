"""Curvilinear normal modes: the modes of TS1 expressed in bonds, angles and dihedrals.

Straight-line (Cartesian) modes read a turning group as a stretch of its bonds, which gives large false displacements to the C-H stretches and
harmonic energies larger than the energy the reaction releases. The same modes taken to first order in internal coordinates do not have that fault.
"""
from __future__ import annotations

import numpy as np

COV = {"H": 0.31, "C": 0.76, "F": 0.57}


def graph(elements, geometries, forming):
    """Bonds of any of the geometries, plus the forming bond, then all angles and dihedrals over those bonds."""
    n = len(elements)
    bonds = {tuple(sorted(forming))}
    for x in geometries:
        d = np.linalg.norm(x[:, None] - x[None], axis=-1)
        for i in range(n):
            for j in range(i):
                if d[i, j] < 1.25 * (COV[elements[i]] + COV[elements[j]]):
                    bonds.add((j, i))
    nbr = {i: [] for i in range(n)}
    for i, j in sorted(bonds):
        nbr[i].append(j)
        nbr[j].append(i)
    angles = [(i, j, k) for j in range(n) for i in nbr[j] for k in nbr[j] if i < k]
    dihed = [(i, j, k, l) for j, k in sorted(bonds) for i in nbr[j] for l in nbr[k] if i != k and l != j and i != l]
    return sorted(bonds), angles, dihed


def values(x, bonds, angles, dihed):
    """Internal coordinates at x (natoms, 3). Lengths in the units of x, angles in radians."""
    r = [np.linalg.norm(x[i] - x[j]) for i, j in bonds]
    a = []
    for i, j, k in angles:
        u, v = x[i] - x[j], x[k] - x[j]
        a.append(np.arctan2(np.linalg.norm(np.cross(u, v)), u @ v))
    d = []
    for i, j, k, l in dihed:
        b0, b1, b2 = x[j] - x[i], x[k] - x[j], x[l] - x[k]
        n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
        d.append(np.arctan2(np.cross(n1, n2) @ b1 / np.linalg.norm(b1), n1 @ n2))
    return np.array(r + a + d)


class Curvilinear:
    """Coordinates Q with xi - xi_TS = (B L) Q to first order at TS1, where B = d xi / d q on mass-weighted q and L are the TS1 modes."""

    def __init__(self, elements, geometries, forming, q0, m, L):
        self.g = graph(elements, geometries, forming)
        self.sm = np.sqrt(m)
        self.xi0 = self.xi(q0)
        self.nb, self.na, self.nd = (len(v) for v in self.g)
        B = np.zeros((len(self.xi0), len(q0)))
        for i in range(len(q0)):
            e = np.zeros(len(q0))
            e[i] = 1e-3
            B[:, i] = self.wrap(self.xi(q0 + e) - self.xi(q0 - e)) / 2e-3
        self.BL = B @ L
        self.pinv = np.linalg.pinv(self.BL, rcond=1e-6)
        self.cond = np.linalg.cond(self.BL)

    def xi(self, q):
        return values((q / self.sm).reshape(-1, 3), *self.g)

    def wrap(self, d):
        d = d.copy()
        if self.nd:
            d[-self.nd:] = (d[-self.nd:] + np.pi) % (2 * np.pi) - np.pi
        return d

    def __call__(self, q):
        """Curvilinear mode coordinates of the geometry q, relative to TS1. The result depends on q only through its internal coordinates."""
        return self.pinv @ self.wrap(self.xi(q) - self.xi0)
