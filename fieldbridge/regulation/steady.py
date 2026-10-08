"""Steady states: reached by integrating the drift, refined by Newton's method, with their stability.

Conservation laws L (rows with L F = 0) fix the totals L q of the initial state; the Newton step then solves the
stacked system [J; L] dq = -[F; L q - totals], and stability is judged on the Jacobian restricted to the null space
of L (the directions the dynamics can move along).
"""
from __future__ import annotations

from types import SimpleNamespace
from typing import Callable, Dict, Optional

import numpy as np
from scipy.integrate import LSODA, solve_ivp
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


FLOOR = 1e-6           # Newton step / |q| within which a residual that has stopped falling is the integrator's error


def at_floor(m: Regulated, p: np.ndarray, q: np.ndarray, L: Optional[np.ndarray], res: float,
             res_before: Optional[float]) -> bool:
    """True when a trajectory has reached the error floor of the integrator: its residual res did not fall over the
    last chunk (from res_before), at a stable state whose Newton step to the steady state is below FLOOR of |q|.
    A residual set by the dynamics keeps falling, by e^-10 or more over a chunk of ten slowest relaxation times and
    by the growth of the chunk even for an algebraic approach (an integrator winding up without a steady state);
    the error of the integrator does not. A settle test below that error (|F| <= 1e-9 fast |q| against rtol 1e-10)
    passes or fails with the rounding of the platform, and a run that fails integrates all its chunks, each 1.5
    times longer."""
    if res_before is None or res < res_before:
        return False
    lam = rates(m, p, q, L)
    if lam.size and not lam.real.max() < -1e-10 * float(np.abs(lam).max()):
        return False
    J, F = m.jac(q, p), m.f(q, p)
    A, b = (np.vstack([J, L]), np.concatenate([F, np.zeros(len(L))])) if L is not None and len(L) else (J, F)
    dq = np.linalg.lstsq(A, b, rcond=None)[0]
    return bool(np.linalg.norm(dq) <= FLOOR * max(np.linalg.norm(q), 1e-300))


BUDGET = 100_000       # right-hand-side calls of LSODA in one chunk, after which BDF integrates the chunk


def _lsoda(fun: Callable, jac: Optional[Callable], span, x0: np.ndarray, rtol: float, atol: float,
           t_eval: np.ndarray) -> Optional[SimpleNamespace]:
    """LSODA stepped as solve_ivp steps it, with the same values at t_eval; None once it has used more than BUDGET
    right-hand-side calls. (An exception raised in the right-hand side would stop it too, but with the Fortran LSODA
    of older SciPy, 1.14 here, it aborts the interpreter.)"""
    solver = LSODA(fun, float(span[0]), x0, float(span[1]), rtol=rtol, atol=atol, jac=jac)
    ts, ys, i, message = [], [], 0, None
    while solver.status == "running":
        message = solver.step()
        if solver.status == "failed":
            break
        if solver.nfev > BUDGET:
            return None
        j = int(np.searchsorted(t_eval, solver.t, side="right"))
        if j > i:
            ts.append(t_eval[i:j])
            ys.append(solver.dense_output()(t_eval[i:j]))
            i = j
    return SimpleNamespace(t=np.hstack(ts) if ts else np.zeros(0), y=np.hstack(ys) if ys else np.zeros((len(x0), 0)),
                           success=solver.status != "failed", message=message, nfev=solver.nfev)


def integrate_chunk(fun: Callable, jac: Callable, span, x0: np.ndarray, rtol: float, atol: float, t_eval: np.ndarray,
                    lsoda_jac: bool = True, stiff: bool = False):
    """One chunk: LSODA, or BDF with the analytic Jacobian jac when stiff is set or LSODA needs more than BUDGET
    right-hand-side calls. Returns the solution and whether BDF integrated the chunk.

    LSODA chooses between its non-stiff and stiff methods by estimates that fail on a response that keeps getting
    stiffer (an integrator winding up at the edge of its range, qian2018_quasi_control2 at d = 10): the same chunk
    then takes from 2e4 to more than 3e6 calls with a change of the rounding. BDF takes 6e3 there, but it is slower on
    most responses and crawls at a kink of the drift (the abs terms of chemotaxis_tu2008), where LSODA does not."""
    if not stiff:
        sol = _lsoda(fun, jac if lsoda_jac else None, span, x0, rtol, atol, t_eval)
        if sol is not None:
            return sol, False
    return solve_ivp(fun, span, x0, method="BDF", jac=jac, rtol=rtol, atol=atol, t_eval=t_eval), True


def output_along(m: Regulated, p: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """The output y at the states Q (one per row), in one vectorized call."""
    return m.y_batch(Q, np.repeat(p[None, :], len(Q), axis=0))


def settle(m: Regulated, p: np.ndarray, q0: np.ndarray, L: Optional[np.ndarray] = None, tol: float = 1e-9,
           max_chunks: int = 14) -> Dict[str, object]:
    """Integrate from q0 in chunks until |F| is small, or the residual has stopped falling at the integrator's error
    floor (Newton's method then refines the state); report the oscillation of y over the last chunk. The total time
    is bounded (chunks of 10, 15, 22, ... slowest relaxation times); a trajectory that grows by more than 1e8 is
    reported as diverging, and one whose output keeps oscillating with an undiminished amplitude is stopped and
    reported as oscillating."""
    q = np.array(q0, float)
    size0 = max(float(np.linalg.norm(q)), 1.0)
    t_total, chunk = 0.0, _time_scale(m, p, q, L)
    jac = (lambda t, x: m.jac(x, p))
    last, previous, res_before, stiff = None, None, None, False
    for _ in range(max_chunks):
        ts = np.linspace(0.0, chunk, 400)
        sol, stiff = integrate_chunk(lambda t, x: m.f(x, p), jac, (0.0, chunk), q, 1e-10,
                                     1e-12 * max(1.0, float(np.max(np.abs(q)))), ts, stiff=stiff)
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
        res = float(np.linalg.norm(F)) / max(fast * max(np.linalg.norm(q), 1e-12), 1e-300)
        if np.linalg.norm(F) <= tol * fast * max(np.linalg.norm(q), 1e-12) or at_floor(m, p, q, L, res, res_before):
            return {"q": q, "converged": True, "time": t_total, "oscillation": last}
        res_before = res
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
