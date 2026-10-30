"""Friction from explicit modes: the minimal surface coupled to harmonic bath modes instead of a Langevin bath.

Each of the two coordinates is coupled to its own set of M harmonic modes (Caldeira-Leggett form, with the
counterterm, so that the potential of mean force is the surface itself):

    V_bath = sum_j (1/2) w_j^2 (Q_j - c_j x / w_j^2)^2 + sum_j (1/2) w_j^2 (R_j - c_j y / w_j^2)^2,

    w_j = j * w_max / M,   c_j^2 = (2/pi) gamma w_j^2 (w_max / M),

a discretized Ohmic spectral density J(w) = gamma w up to w_max. With many modes the memory kernel
K(t) = sum_j c_j^2 / w_j^2 cos(w_j t) tends to 2 gamma delta(t) and the dynamics to Langevin friction gamma; with few
modes it recurs at t = 2 pi M / w_max and the modes exchange energy coherently. The whole system is conservative: the
only noise is the thermal distribution of the modes at the start.

The law for this system is exact in the linear stage, because the modes are linear: the transverse coordinate and its
M modes form a linear system along the reference path x_ref(t), which is itself integrated with the x modes (at rest
at the start). The decision is read at the largest transverse variance beyond the VRI during the first excursion, as in
ptsb.linear_law.
"""
from __future__ import annotations

from math import pi, sqrt
from typing import Dict, List

import numpy as np
from scipy.integrate import solve_ivp

from ptsb import Phi, Minimal, ts_frame


def modes(gamma: float, M: int, w_max: float = 4.0):
    w = np.arange(1, M + 1) * w_max / M
    c = np.sqrt(2 / pi * gamma * w ** 2 * (w_max / M))
    return w, c


def sample(s, n: int, rng, kT: float, gamma: float, M: int, w_max: float = 4.0, speed=None):
    """Thermal start on the dividing line through TS1 (as ptsb.sample_ts) with the modes in equilibrium about the
    displaced coordinate: Q_j ~ N(c_j x / w_j^2, kT / w_j^2), P_j ~ N(0, kT)."""
    from ptsb import sample_ts
    x, y, px, py = sample_ts(s, n, rng, kT=kT, speed=speed)
    w, c = modes(gamma, M, w_max)
    Q = (c / w ** 2)[None, :] * x[:, None] + rng.normal(0, 1, (n, M)) * np.sqrt(kT) / w[None, :]
    R = (c / w ** 2)[None, :] * y[:, None] + rng.normal(0, 1, (n, M)) * np.sqrt(kT) / w[None, :]
    PQ = rng.normal(0, sqrt(kT), (n, M))
    PR = rng.normal(0, sqrt(kT), (n, M))
    return x, y, px, py, Q, R, PQ, PR


def _ops(s, backend):
    """Force and assignment of the surface in either backend: (grad(x, y), check(x, y) -> inside, side, back)."""
    if isinstance(s, Minimal):
        from engine import coeffs, geometry, grad_poly, V_poly
        A = [float(a) for a in coeffs(s)]
        level, xc, qx, qy, nx, ny, xback = (float(v) for v in geometry(s))

        def check(x, y):
            return (x > xc) & (V_poly(A, x, y) < level), (x - qx) * nx + (y - qy) * ny, \
                (x < xback) | (abs(x) > 5.0) | (abs(y) > 5.0)
        return (lambda x, y: grad_poly(A, x, y)), check
    # CurvedValley, kappa > 0 (no pseudo-minimum terms in these tests)
    if backend == "torch":
        import torch
        sqrt_, atan2 = torch.sqrt, torch.atan2
    else:
        sqrt_, atan2 = np.sqrt, np.arctan2
    R, d, eps = 1.0 / s.kappa, s.d, s.eps
    level, sc = float(s.commit_level()), float(s.s_commit())
    nts2, sback = float(s.points_sn()["TS2"][1]), float(s.points_sn()["TS1"][0] - 0.4)

    def grad(x, y):
        rho2 = x * x + (R - y) ** 2
        rho = sqrt_(rho2)
        ss, n = R * atan2(x, R - y), R - rho
        Vs = 8 * ss * ss - 8 * ss + n * n * (n * n - 2)
        Vn = n + (ss - d) * (4 * n ** 3 - 4 * n) + eps
        return Vs * R * (R - y) / rho2 - Vn * x / rho, Vs * R * x / rho2 + Vn * (R - y) / rho

    def check(x, y):
        ss, n = R * atan2(x, R - y), R - sqrt_(x * x + (R - y) ** 2)
        V = 8 / 3 * ss ** 3 - 4 * ss * ss + 0.5 * n * n + (ss - d) * n * n * (n * n - 2) + eps * n
        return (ss > sc) & (V < level), n - nts2, (ss < sback) | (abs(n) > 3.0) | (ss > 4.0)
    return grad, check


def run(s, init, gamma: float, M: int, w_max: float = 4.0, dt: float = 0.005, t_max: float = 60.0,
        check_every: int = 4, backend: str = "numpy", device: str = "cuda") -> np.ndarray:
    """Velocity Verlet for the surface and its modes; first-entry assignment as in ptsb.run (+1 / -1 product side,
    0 returned behind TS1 or escaped, 2 undecided). M = 0 and gamma = 0 runs the bare surface."""
    w, c = modes(gamma, M, w_max)
    x, y, px, py, Q, R, PQ, PR = init
    if backend == "torch":
        import torch
        dev = torch.device(device)
        T = lambda a: torch.as_tensor(np.asarray(a, float), dtype=torch.float64, device=dev)
        x, y, px, py, Q, R, PQ, PR, w, c = (T(v) for v in (x, y, px, py, Q, R, PQ, PR, w, c))
        full = lambda v: torch.full(x.shape, v, dtype=torch.int64, device=dev)
        where, sign = torch.where, torch.sign
        to_np = lambda v: v.cpu().numpy()
        as_long = lambda v: v.long()
    else:
        x, y, px, py = (np.array(v, float) for v in (x, y, px, py))
        Q, R, PQ, PR = (np.array(v, float) for v in (Q, R, PQ, PR))
        full = lambda v: np.full(x.shape, v, np.int64)
        where, sign = np.where, np.sign
        to_np = lambda v: v
        as_long = lambda v: v.astype(np.int64)
    grad, check = _ops(s, backend)
    w2 = w ** 2
    cw = c / w2

    def forces(x, y, Q, R):
        gx, gy = grad(x, y)
        dQ = Q - cw * x[:, None]
        dR = R - cw * y[:, None]
        return -gx + (c * dQ).sum(1), -gy + (c * dR).sum(1), -w2 * dQ, -w2 * dR

    out = full(2)
    alive = full(1) > 0
    fx, fy, fQ, fR = forces(x, y, Q, R)
    for k in range(int(round(t_max / dt))):
        px = px + 0.5 * dt * fx
        py = py + 0.5 * dt * fy
        PQ = PQ + 0.5 * dt * fQ
        PR = PR + 0.5 * dt * fR
        x = x + dt * px
        y = y + dt * py
        Q = Q + dt * PQ
        R = R + dt * PR
        fx, fy, fQ, fR = forces(x, y, Q, R)
        px = px + 0.5 * dt * fx
        py = py + 0.5 * dt * fy
        PQ = PQ + 0.5 * dt * fQ
        PR = PR + 0.5 * dt * fR
        if (k + 1) % check_every == 0:
            ins, side, bk = check(x, y)
            inside = alive & ins
            sd = as_long(sign(side))
            back = alive & bk
            out = where(back, full(0), out)
            out = where(inside, sd, out)
            alive = alive & ~back & ~inside
            if not bool(alive.any()):
                break
            m = alive.double() if backend == "torch" else alive.astype(float)
            px, py, fx, fy = px * m, py * m, fx * m, fy * m
            PQ, PR, fQ, fR = PQ * m[:, None], PR * m[:, None], fQ * m[:, None], fR * m[:, None]
    return to_np(out)


def law(s: Minimal, speed: float, kT: float, gamma: float, M: int, w_max: float = 4.0, t_end: float = 14.0,
        n_eval: int = 2801) -> Dict[str, float]:
    """The exact linear-stage law with the modes (no Markov approximation)."""
    w, c = modes(gamma, M, w_max)
    w2 = w ** 2
    K0 = float((c ** 2 / w2).sum())
    _, _, _, om, _ = ts_frame(s)
    D = 2 + 2 * M
    vy0 = kT / om ** 2
    y0 = -s.grad(0.0, 0.0)[1] / s.hessian(0.0, 0.0)[1, 1]
    # initial mean and covariance of (y, p, R_1..R_M, PR_1..PR_M)
    m0 = np.r_[y0, 0.0, c / w2 * y0, np.zeros(M)]
    S0 = np.zeros((D, D))
    g = np.r_[1.0, 0.0, c / w2, np.zeros(M)]           # R_j = (c_j / w_j^2) y + xi_j
    S0 += vy0 * np.outer(g, g)
    S0[1, 1] = kT
    S0[2:2 + M, 2:2 + M] += np.diag(kT / w2)
    S0[2 + M:, 2 + M:] = np.diag(np.full(M, kT))

    def Amat(xr):
        A = np.zeros((D, D))
        A[0, 1] = 1
        A[1, 0] = -(s.hessian(xr, 0.0)[1, 1] + K0)
        A[1, 2:2 + M] = c
        A[2:2 + M, 2 + M:] = np.eye(M)
        A[2 + M:, 0] = c
        A[2 + M:, 2:2 + M] = -np.diag(w2)
        return A

    def rhs(t, z):
        xr, v = z[0], z[1]
        Qm, Pm = z[2:2 + M], z[2 + M:2 + 2 * M]
        o = 2 + 2 * M
        m = z[o:o + D]
        S = z[o + D:].reshape(D, D)
        A = Amat(xr)
        b = np.zeros(D)
        b[1] = -s.grad(xr, 0.0)[1]
        acc = -s.grad(xr, 0.0)[0] + (c * (Qm - c / w2 * xr)).sum()
        return np.r_[v, acc, Pm, -w2 * Qm + c * xr, A @ m + b, (A @ S + S @ A.T).ravel()]

    def behind(t, z):
        return z[0] + 0.4
    behind.terminal = True

    z0 = np.r_[0.0, speed, np.zeros(2 * M), m0, S0.ravel()]
    t = np.linspace(0, t_end, n_eval)
    sol = solve_ivp(rhs, (0, t_end), z0, t_eval=t, rtol=1e-6, atol=1e-9, method="DOP853", events=(behind,))
    o = 2 + 2 * M
    xr = sol.y[0]
    my = sol.y[o]
    vy = sol.y[o + D]                                    # S[0, 0]
    beyond = xr > s.vri()
    if not beyond.any():
        return {"P_law": float("nan")}
    k0 = int(np.argmax(beyond))
    back = np.where(~beyond[k0:])[0]
    k1 = k0 + (int(back[0]) if back.size else len(xr) - k0)
    kb = k0 + int(np.argmax(vy[k0:k1]))
    z = my[kb] / sqrt(vy[kb])
    return {"P_law": float(Phi(z)), "z": float(z), "x_decision": float(xr[kb]),
            "growth": float(sqrt(vy[kb] / vy0))}


def law_thermal(s: Minimal, kT: float, gamma: float, M: int, n_speeds: int = 16, **kw) -> Dict[str, float]:
    u = (np.arange(n_speeds) + 0.5) / n_speeds
    v = np.sqrt(-2.0 * kT * np.log(1.0 - u))
    L = [law(s, float(vi), kT, gamma, M, **kw) for vi in v]
    return {"P_law": float(np.nanmean([l["P_law"] for l in L])),
            "growth": float(np.nanmean([l.get("growth", np.nan) for l in L]))}
