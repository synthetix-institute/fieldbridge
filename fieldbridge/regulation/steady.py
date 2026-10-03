"""Steady states: reached by integrating the drift, refined by Newton's method, with their stability.

Conservation laws L (rows with L F = 0) fix the totals L q of the initial state; the Newton step then solves the
stacked system [J; L] dq = -[F; L q - totals], and stability is judged on the Jacobian restricted to the null space
of L (the directions the dynamics can move along).
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import null_space

from .spec import Regulated


def rates(m: Regulated, p: np.ndarray, q: np.ndarray, L: Optional[np.ndarray] = None) -> np.ndarray:
    """Eigenvalues of the Jacobian, restricted to the directions allowed by the conservation laws."""
    J = m.jac(q, p)
    if L is not None and len(L):
        U = null_space(L)
        J = U.T @ J @ U
    return np.linalg.eigvals(J) if J.size else np.zeros(0)


def _time_scale(m: Regulated, p: np.ndarray, q: np.ndarray, L) -> float:
    lam = rates(m, p, q, L)
    re = np.abs(lam.real)
    re = re[re > 1e-9 * max(re.max(initial=0.0), 1e-300)] if re.size else re
    slow = re.min() if re.size else 1.0
    return 10.0 / max(slow, 1e-12)


def output_along(m: Regulated, p: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """The output y at the states Q (one per row), in one vectorized call."""
    return m.y_batch(Q, np.repeat(p[None, :], len(Q), axis=0))


def settle(m: Regulated, p: np.ndarray, q0: np.ndarray, L: Optional[np.ndarray] = None, tol: float = 1e-9,
           max_chunks: int = 14) -> Dict[str, object]:
    """Integrate from q0 in chunks until |F| is small (Newton's method then refines the state); report the
    oscillation of y over the last chunk. The total time is bounded (chunks of 10, 15, 22, ... slowest relaxation
    times); a trajectory that grows by more than 1e8 is reported as diverging, and one whose output keeps
    oscillating with an undiminished amplitude is stopped and reported as oscillating."""
    q = np.array(q0, float)
    size0 = max(float(np.linalg.norm(q)), 1.0)
    t_total, chunk = 0.0, _time_scale(m, p, q, L)
    jac = (lambda t, x: m.jac(x, p))
    last, previous = None, None
    for _ in range(max_chunks):
        ts = np.linspace(0.0, chunk, 400)
        sol = solve_ivp(lambda t, x: m.f(x, p), (0.0, chunk), q, method="LSODA", jac=jac, rtol=1e-10,
                        atol=1e-12 * max(1.0, float(np.max(np.abs(q)))), t_eval=ts)
        if not sol.success or not np.all(np.isfinite(sol.y[:, -1])):
            return {"q": q, "converged": False, "time": t_total, "failure": sol.message, "oscillation": None}
        t_total += chunk
        ys = output_along(m, p, sol.y.T)
        previous, last = last, float(ys.max() - ys.min())
        q = sol.y[:, -1]
        if np.linalg.norm(q) > 1e8 * size0:
            return {"q": q, "converged": False, "time": t_total, "failure": "diverging", "oscillation": last}
        F = m.f(q, p)
        lam = np.abs(rates(m, p, q, L).real)
        fast = float(lam.max()) if lam.size else 1.0
        if np.linalg.norm(F) <= tol * fast * max(np.linalg.norm(q), 1e-12):
            return {"q": q, "converged": True, "time": t_total, "oscillation": last}
        scale_y = max(abs(float(ys[-1])), 1e-12)
        if previous is not None and last > 1e-6 * scale_y and last > 0.9 * previous:
            return {"q": q, "converged": False, "time": t_total, "failure": "oscillating", "oscillation": last}
        chunk *= 1.5
    return {"q": q, "converged": False, "time": t_total, "oscillation": last}


def refine(m: Regulated, p: np.ndarray, q: np.ndarray, L: Optional[np.ndarray] = None,
           totals: Optional[np.ndarray] = None, iters: int = 60) -> Dict[str, object]:
    """Newton's method on F = 0 (stacked with L q = totals when there are conservation laws)."""
    q = np.array(q, float)
    for it in range(iters):
        F, J = m.f(q, p), m.jac(q, p)
        if L is not None and len(L):
            A, b = np.vstack([J, L]), np.concatenate([F, L @ q - totals])
        else:
            A, b = J, F
        dq = np.linalg.lstsq(A, b, rcond=None)[0]
        q = q - dq
        if np.linalg.norm(dq) <= 1e-14 * max(np.linalg.norm(q), 1e-300):
            break
    F = m.f(q, p)
    return {"q": q, "residual": float(np.linalg.norm(F)), "iterations": it + 1}


def steady_state(m: Regulated, p: np.ndarray, q0: np.ndarray, L: Optional[np.ndarray] = None) -> Dict[str, object]:
    """The steady state reached from q0: converged, stable, the state, its output and the rates around it."""
    s = settle(m, p, q0, L)
    totals = L @ np.asarray(q0, float) if L is not None and len(L) else None
    r = refine(m, p, s["q"], L, totals)
    moved = float(np.linalg.norm(r["q"] - s["q"]) / max(np.linalg.norm(s["q"]), 1e-300))
    lam0 = rates(m, p, r["q"], L)
    fast = float(np.abs(lam0).max()) if lam0.size else 1.0
    equilibrium = bool(np.all(np.isfinite(r["q"])) and r["residual"] < 1e-8 * fast * max(1.0, np.linalg.norm(r["q"])))
    # keep Newton's point unless it jumped away from a trajectory that had settled
    q = r["q"] if equilibrium and (moved < 1e-4 or not s["converged"]) else s["q"]
    lam = rates(m, p, q, L)
    fast = float(np.abs(lam).max()) if lam.size else 1.0
    stable = bool(lam.size == 0 or lam.real.max() < -1e-10 * fast)
    return {"q": q, "y": m.y(q, p), "converged": bool(s["converged"]), "equilibrium": equilibrium,
            "stable": stable, "rates": lam, "time": s["time"], "oscillation": s.get("oscillation"),
            "newton_moved": moved}
