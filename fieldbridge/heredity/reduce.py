"""The reduction of a body at its threshold along the size (letter K).

The symmetric state q_sym of the body loses stability at a pitchfork at the size L_c: the critical eigenvalue
lambda(L) of the Jacobian at q_sym crosses zero. With the right and left eigenvectors v and w (w.v = 1) the critical
coordinate is s = w.(q - q_sym) and

    ds/dt = a (L - L_c) s - b s^3 + sqrt(2 D_s) xi,   a = d lambda/dL at L_c,   D_s = w^T D w,

where b comes from the drift along the slow manifold (the other directions slaved to s at L_c): b > 0 is a supercritical
pitchfork, for which the law of inheritance holds.
"""
from __future__ import annotations

from typing import Dict

import numpy as np
from scipy.optimize import brentq, root


def reduce_at_threshold(body, q_sym0: np.ndarray, L_lo: float, L_hi: float, s_max: float = 0.02) -> Dict:
    """L_c of the symmetric state (continued from q_sym0) between L_lo and L_hi, with a, b, v and w."""
    cache: Dict[float, np.ndarray] = {}

    def sym(L):
        key = round(L, 12)
        if key not in cache:
            near = min(cache.items(), key=lambda kv: abs(kv[0] - L))[1] if cache else q_sym0
            cache[key] = body.steady(near, L)
        return cache[key]

    def lam(L):
        ev = np.linalg.eigvals(body.jac(sym(L), L))
        return float(ev[np.argmax(ev.real)].real)

    lo, hi = lam(L_lo), lam(L_hi)
    if not (lo < 0 < hi):
        raise ValueError(f"no threshold between {L_lo} and {L_hi}: the largest eigenvalue goes from {lo:.3g} to "
                         f"{hi:.3g}")
    Lc = brentq(lam, L_lo, L_hi, xtol=1e-12)
    q = sym(Lc)
    J = body.jac(q, Lc)
    ev, V = np.linalg.eig(J)
    v = np.real(V[:, int(np.argmin(np.abs(ev)))])
    v /= np.linalg.norm(v)
    evl, W = np.linalg.eig(J.T)
    w = np.real(W[:, int(np.argmin(np.abs(evl)))])
    w /= w @ v
    dL = 1e-5 * max(1.0, abs(Lc))
    a = (lam(Lc + dL) - lam(Lc - dL)) / (2 * dL)
    n = body.n
    P = np.eye(n) - np.outer(v, w)
    U = np.linalg.svd(P)[0][:, : n - 1] if n > 1 else np.zeros((1, 0))
    ss = np.linspace(-s_max, s_max, 9)
    ss = ss[ss != 0]
    fs = []
    for s in ss:
        if n > 1:
            sol = root(lambda y: U.T @ (P @ body.f(q + s * v + U @ y, Lc)[0]), np.zeros(n - 1), tol=1e-13)
            qq = q + s * v + U @ sol.x
        else:
            qq = q + s * v
        fs.append(float(w @ body.f(qq, Lc)[0]))
    coef = np.linalg.lstsq(np.vstack([ss, ss ** 3, ss ** 5]).T, np.array(fs), rcond=None)[0]
    return {"L_c": float(Lc), "q_sym": q, "a": float(a), "b": float(-coef[1]), "v": v, "w": w,
            "linear_at_threshold": float(coef[0])}


def critical_eigenvalue(body, red: Dict, L: float) -> float:
    """The eigenvalue of the critical mode at the symmetric state of size L: the one whose eigenvector overlaps v the
    most (the symmetric state is followed from L_c)."""
    q = body.steady(red["q_sym"], L) if hasattr(body, "steady") else red["q_sym"]
    ev, V = np.linalg.eig(body.jac(q, L))
    Vr = np.real(V)
    overlap = np.abs(Vr.T @ red["v"]) / np.linalg.norm(Vr, axis=0).clip(1e-300)
    return float(ev[int(np.argmax(overlap))].real)
