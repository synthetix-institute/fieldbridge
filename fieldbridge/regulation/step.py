"""Response of the output to a step of the input, and the calibration of an integrator.

The state is integrated from the steady state at the reference input together with two running integrals,

    I1 = int (y - y_ref) dt        and, for an integrator phi,   I2 = int g(q) (y - y0) dt,

so that both come with the accuracy of the solver. For an integrator the identity I2 = phi(end) - phi(start) holds
exactly (for stage 3, I2 = int w.F dt); the ratio I2 / delta phi is a calibration of the simulation and of the
integrator found, not a law.
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np
from scipy.integrate import solve_ivp

from .integrator import phi, rate
from .spec import Regulated
from .steady import output_along, rates, refine


def step_response(m: Regulated, p0: np.ndarray, q0: np.ndarray, u1: float, integ: Optional[Dict] = None,
                  L: Optional[np.ndarray] = None, tol: float = 1e-9, max_chunks: int = 14) -> Dict[str, object]:
    p1 = p0.copy()
    p1[m.u_index] = u1
    y_ref = m.y(q0, p0)
    has_int = bool(integ and integ.get("found"))
    y0 = float(integ["set_point"]) if has_int else y_ref

    def rhs(t, x):
        q = x[: m.n]
        y = m.y(q, p1)
        extra = [y - y_ref]
        if has_int:
            extra.append(rate(integ, m, q, p1, y))
        return np.concatenate([m.f(q, p1), extra])

    x = np.concatenate([np.asarray(q0, float), np.zeros(2 if has_int else 1)])
    lam = np.abs(rates(m, p1, q0, L).real)
    lam = lam[lam > 1e-9 * max(lam.max(initial=0.0), 1e-300)]
    chunk = 10.0 / max(lam.min() if lam.size else 1.0, 1e-12)
    t0, ts_all, ys_all, settled = 0.0, [], [], False
    for _ in range(max_chunks):
        ts = np.linspace(t0, t0 + chunk, 2000)
        sol = solve_ivp(rhs, (t0, t0 + chunk), x, method="LSODA", rtol=1e-11,
                        atol=1e-13 * max(1.0, float(np.max(np.abs(x[: m.n])))), t_eval=ts)
        if not sol.success:
            return {"u1": u1, "failure": sol.message}
        ts_all.append(ts)
        ys_all.append(output_along(m, p1, sol.y[: m.n].T))
        x, t0 = sol.y[:, -1], t0 + chunk
        q = x[: m.n]
        F = m.f(q, p1)
        fast = float(np.abs(rates(m, p1, q, L)).max(initial=1.0))
        if np.linalg.norm(F) <= tol * fast * max(np.linalg.norm(q), 1e-12):
            settled = True
            break
        chunk *= 1.5
    t, y = np.concatenate(ts_all), np.concatenate(ys_all)
    totals = L @ np.asarray(q0, float) if L is not None and len(L) else None
    fine = refine(m, p1, x[: m.n], L, totals)
    q_end = fine["q"] if np.linalg.norm(fine["q"] - x[: m.n]) < 1e-6 * max(np.linalg.norm(x[: m.n]), 1e-12) \
        else x[: m.n]
    y_end = m.y(q_end, p1)
    dev = y - y_ref
    i_peak = int(np.argmax(np.abs(dev)))
    peak = float(dev[i_peak])
    final = float(y_end - y_ref)
    # return time: last time the output is farther than 1% of the peak deviation from its final value
    far = np.nonzero(np.abs(y - y_end) > 0.01 * abs(peak))[0] if peak != 0 else np.array([], int)
    out = {"u1": float(u1), "settled": settled, "y_ref": float(y_ref), "y_end": float(y_end), "peak": peak,
           "t_peak": float(t[i_peak]), "final": final,
           "final_over_peak": final / peak if peak != 0 else None,
           "return_time": float(t[far[-1]]) if far.size else 0.0, "integral_of_deviation": float(x[m.n]),
           "q_end": x[: m.n].tolist()}
    if has_int:
        dphi = phi(integ["coefficients"], m.variables, x[: m.n]) - phi(integ["coefficients"], m.variables, q0)
        out["delta_phi"] = float(dphi)
        out["integral_of_gain_times_error"] = float(x[m.n + 1])
        out["calibration_ratio"] = float(x[m.n + 1] / dphi) if abs(dphi) > 1e-300 else None
    return out
