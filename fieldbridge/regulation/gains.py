"""Static gains of the output, with and without a clamped variable, and the role of each variable.

At a steady state F(q*, u) = 0 the output y = h(q*, u) changes with the input as

    G = dy/du = h_u - grad h . J^-1 F_u = det[[J, F_u], [grad h, h_u]] / det J,

the second form being the Schur complement of the bordered Jacobian. With conservation laws L the shift of the
steady state solves the stacked system [J; L] dq = [-F_u; 0]. The relative sensitivity S = (u/y) G is the
dimensionless gain; perfect adaptation is S = 0.

Clamping a variable (slot C) holds it at its steady value, as a reservoir or a chemostat would, and removes its
equation: the open-loop gain through that variable. The ratio rho = G_closed / G_open is the fraction of the step
that remains; for proportional feedback with loop gain Lg it is 1/(1 + Lg).
"""
from __future__ import annotations

from collections import deque
from typing import Dict, List, Optional, Set

import numpy as np

from .spec import Regulated


def static_gain(m: Regulated, p: np.ndarray, q: np.ndarray, L: Optional[np.ndarray] = None) -> Dict[str, object]:
    J, Fu, gh, hu = m.jac(q, p), m.fu(q, p), m.gradh(q, p), m.hu(q, p)
    bordered = None
    if L is not None and len(L):
        dq = np.linalg.lstsq(np.vstack([J, L]), np.concatenate([-Fu, np.zeros(len(L))]), rcond=None)[0]
    else:
        dq = -np.linalg.solve(J, Fu)
        detJ = np.linalg.det(J)
        if detJ != 0:
            bordered = float(np.linalg.det(np.block([[J, Fu[:, None]], [gh[None, :], np.array([[hu]])]])) / detJ)
    G = float(hu + gh @ dq)
    y, u = m.y(q, p), float(p[m.u_index])
    # S is reported only where y is not zero within the size of the terms that make it up
    size = abs(hu * u) + float(np.abs(gh * q).sum())
    S = G * u / y if u != 0 and abs(y) > 1e-9 * max(size, 1e-300) else None
    return {"G": G, "S": S, "dq_du": dq, "bordered": bordered}


def reference_gain(m: Regulated, p: np.ndarray, q: np.ndarray, L: Optional[np.ndarray] = None) -> float:
    """The scale against which a vanishing gain is judged: the largest open-loop gain over the clamps of each variable
    the output does not read and of all of them together (the local response of the output); without such a
    variable, the sum of the absolute contributions to G."""
    others = [k for k in range(m.n) if not m.output_reads[k]]
    refs = []
    for k in others + ([set(others)] if len(others) > 1 else []):
        c = clamped_gain(m, p, q, k, L)
        if c["G"] is not None:
            refs.append(abs(c["G"]))
    if refs and max(refs) > 0:
        return max(refs)
    g = static_gain(m, p, q, L)
    return abs(m.hu(q, p)) + float(np.abs(m.gradh(q, p) * g["dq_du"]).sum())


def clamped_gain(m: Regulated, p: np.ndarray, q: np.ndarray, k, L: Optional[np.ndarray] = None) -> Dict[str, object]:
    """Gain with variable k (an index, or a set of indices) held at its value in q. Laws that involve a clamped
    variable no longer hold; the others do."""
    ks = {k} if isinstance(k, (int, np.integer)) else set(k)
    keep = [i for i in range(m.n) if i not in ks]
    J, Fu, gh, hu = m.jac(q, p), m.fu(q, p), m.gradh(q, p), m.hu(q, p)
    Jr, Fur, ghr = J[np.ix_(keep, keep)], Fu[keep], gh[keep]
    Lr = None
    if L is not None and len(L):
        rows = [r for r in L if all(abs(r[j]) < 1e-12 * np.max(np.abs(r)) for j in ks)]
        Lr = np.array([r[keep] for r in rows]) if rows else None
    if Lr is not None and len(Lr):
        dq = np.linalg.lstsq(np.vstack([Jr, Lr]), np.concatenate([-Fur, np.zeros(len(Lr))]), rcond=None)[0]
        from scipy.linalg import null_space
        U = null_space(Lr)
        lam = np.linalg.eigvals(U.T @ Jr @ U) if U.size else np.zeros(0)
    else:
        if Jr.size == 0 or abs(np.linalg.det(Jr)) < 1e-300:
            return {"G": None, "S": None, "stable": False, "singular": True}
        dq = -np.linalg.solve(Jr, Fur)
        lam = np.linalg.eigvals(Jr)
    G = float(hu + ghr @ dq)
    y, u = m.y(q, p), float(p[m.u_index])
    fast = float(np.abs(lam).max()) if lam.size else 1.0
    return {"G": G, "S": G * u / y if y != 0 and u != 0 else None,
            "stable": bool(lam.size == 0 or lam.real.max() < -1e-10 * fast), "singular": False}


# ---------------------------------------------------------------------------------------------------- topology
def _reach(m: Regulated, start: Set[int]) -> Set[int]:
    """Variables reached from `start` along the interaction graph (edge j -> i when F_i contains q_j)."""
    seen, queue = set(start), deque(start)
    while queue:
        j = queue.popleft()
        for i in range(m.n):
            if m.depends[i][j] and i not in seen:
                seen.add(i)
                queue.append(i)
    return seen


def role(m: Regulated, k: int) -> str:
    """'output' for a variable the observable reads, 'feedback' for a variable on a loop through the output,
    'feedforward' for one between the input and the output off any such loop, 'none' otherwise."""
    outs = {j for j in range(m.n) if m.output_reads[j]}
    ins = {i for i in range(m.n) if m.input_enters[i]}
    if k in outs:
        return "output"
    reaches_output = bool(outs & _reach(m, {k}))
    reached_from_output = k in _reach(m, outs)
    if reaches_output and reached_from_output:
        return "feedback"
    if reaches_output and k in _reach(m, ins):
        return "feedforward"
    return "none"


def input_reaches_output(m: Regulated) -> bool:
    ins = {i for i in range(m.n) if m.input_enters[i]}
    outs = {j for j in range(m.n) if m.output_reads[j]}
    return m.output_reads_input or bool(outs & _reach(m, ins))


def attenuation(m: Regulated, p: np.ndarray, q: np.ndarray, L: Optional[np.ndarray] = None) -> List[Dict[str, object]]:
    """For each variable that is not read by the output: its role and the remaining fraction G_closed/G_open."""
    closed = static_gain(m, p, q, L)["G"]
    rows = []
    for k, name in enumerate(m.variables):
        r = role(m, k)
        if r == "output":
            continue
        c = clamped_gain(m, p, q, k, L)
        ratio = closed / c["G"] if c["G"] not in (None, 0.0) else None
        rows.append({"variable": name, "role": r, "G_open": c["G"], "open_loop_stable": c["stable"],
                     "remaining_fraction": ratio})
    return rows
