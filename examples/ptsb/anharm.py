"""Friction from anharmonic coupling: the surface coupled to further modes through cubic terms only.

No mode is coupled linearly to x or y, as at a stationary point of a molecule, where the normal modes are orthogonal.
The modes Q_j (frequencies w_j, parity +1 or -1 under y -> -y) enter through

    V = V0(x, y) + sum_j (1/2) w_j^2 (Q_j - a_j(x, y))^2 + (1/2) x Q^T A Q + (1/2) y Q^T B Q + sum phi Q_i Q_j Q_k,

    a_j(x, y) = (1/2) dxx_j x^2 + dxy_j x y + (1/2) dyy_j y^2.

The first group displaces the mode along the path. Its cubic terms are (1/2) D_j x^2 Q_j, C_j x y Q_j, (1/2) E_j y^2 Q_j
with D = w^2 dxx and so on, and the square keeps the relaxed surface min_Q V at V0: the statics of (x, y) at T = 0 do
not change, only the dynamics do (the reaction-path coupling of Miller, Handy and Adams). These forces are linear in Q
and scale with the displacement along the path, so they survive at T -> 0. The second group makes the mode frequencies
depend on (x, y). Its forces are quadratic in Q and scale with the population of the modes. The sparse phi carries energy
among the modes themselves.

A displaced mode that is fast follows a_j(x(t), y(t)) instead of absorbing energy for good. Its kinetic energy then
gives the surface coordinates a position-dependent inertia, with no dissipation: G_ik = delta_ik + sum_j d_i a_j d_k a_j.
`run_following` integrates this limit (Q_j = a_j exactly). `run` records the energy of following, (1/2) a_j'^2, and the
energy lost by (x, y). The latter reads harmonic modes of known friction correctly (`gamma_from_energy`). The energy in the
frame of the displaced mode, also recorded (dE_left), does not: it counts (1/2) a_j'^2 for a mode that has not yet responded
and reads about twice the known friction, so it is kept as a diagnostic and not used.

Natural size of a coupling. With Q = (hbar / w)^(1/2) q, a cubic constant is the same number in dimensionless
coordinates when it scales as g (w_i w_j w_k)^(1/2), so all tables are g sqrt(w...) zeta with zeta standard normal and
each mode has a few partners (molecular cubic tables are sparse). The strength g is the only parameter. For orientation:
if frequency 1 is 150 cm-1 and kT = 0.03 is room temperature, then hbar = 0.022, the energy unit is about 6900 cm-1 and
a term g sqrt(w) x y Q is a cubic constant of about 13 g cm-1 in dimensionless normal coordinates. An estimate, not a fit.

Predictions, all computed before any trajectory and with nothing fitted:
* `equivalent`: second-order cumulant of the quadratic forces. A pair term a x Q_i Q_j has the force correlation
  a^2 (kT/w_i^2) (kT/w_j^2) cos(w_i t) cos(w_j t), hence by the fluctuation-dissipation relation the kernel
  K(t) = kappa [cos((w_i + w_j) t) + cos((w_i - w_j) t)] with kappa = a^2 kT / (2 w_i^2 w_j^2), and a diagonal term
  (1/2) a x Q_i^2 has K(t) = (a^2 kT / 4 w^4)(1 + cos 2 w t), whose constant part is a static force set by the action of
  the mode. The friction is proportional to kT. The equivalent system has harmonic modes at these frequencies, linear
  coupling and no counterterm, so its potential of mean force softens by K(0), as that of the cubic system does.
* `law`: the linear-stage law for the system above. The modes of C x y Q are linear in (y, Q) along the reference path and
  are kept exactly, the modes of D x^2 Q are integrated with the reference path, E shifts the stiffness through the mean
  of Q - a, and the quadratic forces enter through their equivalent modes.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import pi, sqrt
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from ptsb import Minimal, Phi, sample_ts, ts_frame

FORMS = ("xxQ", "xyQ", "yyQ", "xQQ", "yQQ")
HBAR = 0.022                                      # the orientation above; used only by the quantum occupations


@dataclass
class Tables:
    w: np.ndarray
    dxx: np.ndarray = None
    dxy: np.ndarray = None
    dyy: np.ndarray = None
    A: np.ndarray = None
    B: np.ndarray = None
    tri: List[Tuple[int, int, int, float]] = field(default_factory=list)
    Lx: np.ndarray = None                         # linear coupling without counterterm (the equivalent system)
    Ly: np.ndarray = None
    Kx: float = 0.0                               # added stiffness: the counterterm of Lx, Ly, when one is wanted
    Ky: float = 0.0
    fx: float = 0.0                               # constant mean forces
    fy: float = 0.0
    static_x: Optional[Tuple[np.ndarray, np.ndarray]] = None    # (coefficient, mean action) of the static force
    static_y: Optional[Tuple[np.ndarray, np.ndarray]] = None
    frame: str = "xy"                             # "path": the couplings use (s, n) of the curved valley instead of (x, y)

    def __post_init__(self):
        M = self.w.size
        for k in ("dxx", "dxy", "dyy", "Lx", "Ly"):
            if getattr(self, k) is None:
                setattr(self, k, np.zeros(M))
        for k in ("A", "B"):
            if getattr(self, k) is None:
                setattr(self, k, np.zeros((M, M)))

    @property
    def M(self):
        return self.w.size


def frequencies(M: int, w_min: float = 0.5, w_max: float = 8.0) -> np.ndarray:
    return np.linspace(w_min, w_max, M)


def stable(T: Tables, margin: float = 0.05) -> bool:
    """Every mode stays bound over the region the trajectories visit."""
    W = np.diag(T.w ** 2)
    for x in (-0.4, 0.0, 0.7, 1.4):
        for y in (-1.2, 0.0, 1.2):
            if np.linalg.eigvalsh(W + x * T.A + y * T.B)[0] < margin * T.w[0] ** 2:
                return False
    return True


def draw(forms, g: float, M: int, seed: int, partners: int = 3, tries: int = 200, w_min: float = 0.5,
         w_max: float = 8.0) -> Tuple[Optional[Tables], int]:
    """A random sparse table of the requested cubic forms, respecting y -> -y (even modes: odd index). A table in which a
    mode loses its binding along the path is redrawn; returns (table or None, redraws)."""
    w = frequencies(M, w_min, w_max)
    even = np.arange(M) % 2 == 1
    rng = np.random.default_rng([seed, M, int(round(g * 1000))] + [FORMS.index(f) if f in FORMS else 9 for f in forms])
    for attempt in range(tries):
        T = Tables(w=w)
        if "xxQ" in forms:
            T.dxx = np.where(even, g * np.sqrt(w) * rng.normal(0, 1, M), 0.0) / w ** 2
        if "xyQ" in forms:
            T.dxy = np.where(~even, g * np.sqrt(w) * rng.normal(0, 1, M), 0.0) / w ** 2
        if "yyQ" in forms:
            T.dyy = np.where(even, g * np.sqrt(w) * rng.normal(0, 1, M), 0.0) / w ** 2
        if "xQQ" in forms:                         # pairs of the same parity, the diagonal included
            for i in range(M):
                T.A[i, i] += g * w[i] * rng.normal() / 1.0
                same = [j for j in range(M) if j != i and even[j] == even[i]]
                for j in rng.choice(same, size=min(partners, len(same)), replace=False):
                    a = g * sqrt(w[i] * w[j]) * rng.normal() / sqrt(2)      # each pair is drawn from both ends
                    T.A[i, j] += a
                    T.A[j, i] += a
        if "yQQ" in forms:                         # pairs of opposite parity
            for i in range(M):
                other = [j for j in range(M) if even[j] != even[i]]
                for j in rng.choice(other, size=min(partners, len(other)), replace=False):
                    a = g * sqrt(w[i] * w[j]) * rng.normal() / sqrt(2)
                    T.B[i, j] += a
                    T.B[j, i] += a
        if "QQQ" in forms:                         # only among modes with a barrier well above the thermal energy
            hi = [i for i in range(M) if w[i] >= 1.5]
            for _ in range(2 * M):
                i, j, k = rng.choice(hi, 3, replace=False)
                if (even[i] ^ even[j] ^ even[k]):  # an even number of odd modes: a parity-even term
                    T.tri.append((int(i), int(j), int(k), float(g * sqrt(w[i] * w[j] * w[k]) * rng.normal())))
        if stable(T):
            return T, attempt
    return None, tries


def entropic(g: float, M: int, sign: float = 1.0, w_min: float = 0.5, w_max: float = 8.0) -> Tables:
    """Modes that stiffen on one side: (1/2) y sum_j B_jj Q_j^2 with B_jj = sign g w_j, on a surface with no tilt. The mean
    force on y is -(1/2) sum B_jj kT / w_j^2 per unit temperature; each trajectory also carries its own static force,
    since the action of a mode survives the passage."""
    w = frequencies(M, w_min, w_max)
    return Tables(w=w, B=np.diag(sign * g * w))


# ------------------------------------------------------------------------------------- the equivalent harmonic system
def kernel(C: np.ndarray, w: np.ndarray, kT: float):
    """Frequencies and weights of K(t) = sum kappa cos(omega t), the static part of the diagonal apart."""
    om, ka = [], []
    M = w.size
    for i in range(M):
        for j in range(i + 1, M):
            if C[i, j] != 0:
                k = C[i, j] ** 2 * kT / (2 * w[i] ** 2 * w[j] ** 2)
                om += [w[i] + w[j], abs(w[i] - w[j])]
                ka += [k, k]
        if C[i, i] != 0:
            om.append(2 * w[i])
            ka.append(C[i, i] ** 2 * kT / (4 * w[i] ** 4))
    return np.array(om), np.array(ka)


def binned(om: np.ndarray, ka: np.ndarray, bins: int = 16):
    """Merge the kernel into at most `bins` modes spaced evenly in frequency (each at its weighted mean frequency)."""
    if om.size == 0:
        return np.zeros(0), np.zeros(0)
    edges = np.linspace(0.0, om.max() * 1.0001, bins + 1)
    idx = np.digitize(om, edges) - 1
    O, K = [], []
    for b in range(bins):
        m = idx == b
        if m.any():
            O.append(float((om[m] * ka[m]).sum() / ka[m].sum()))
            K.append(float(ka[m].sum()))
    return np.array(O), np.array(K)


def means_and_static(T: Tables, kT: float):
    """Mean forces on x and y from the populations of the modes, and the static forces per trajectory. The static force of
    mode j is -(1/2) a_jj J_j / w_j, with the action J_j exponential of mean kT / w_j."""
    dA, dB = np.diag(T.A), np.diag(T.B)
    fx = -0.5 * float((dA * kT / T.w ** 2).sum())
    fy = -0.5 * float((dB * kT / T.w ** 2).sum())
    return fx, fy


def equivalent(T: Tables, kT: float, bins: int = 32, static: bool = True) -> Tables:
    """The system with the quadratic forces replaced by harmonic modes of the same K(t). Displacement terms are kept."""
    ox, kx = binned(*kernel(T.A, T.w, kT), bins)
    oy, ky = binned(*kernel(T.B, T.w, kT), bins)
    disp = [i for i in range(T.M) if T.dxx[i] or T.dxy[i] or T.dyy[i]]
    w = np.r_[T.w[disp], ox, oy]
    nd, nx, ny = len(disp), ox.size, oy.size
    Lx, Ly = np.zeros(w.size), np.zeros(w.size)
    Lx[nd:nd + nx] = ox * np.sqrt(kx)              # c^2 / omega^2 = kappa
    Ly[nd + nx:] = oy * np.sqrt(ky)
    fx, fy = means_and_static(T, kT)
    S = Tables(w=w, dxx=np.r_[T.dxx[disp], np.zeros(nx + ny)], dxy=np.r_[T.dxy[disp], np.zeros(nx + ny)],
               dyy=np.r_[T.dyy[disp], np.zeros(nx + ny)], Lx=Lx, Ly=Ly, fx=fx, fy=fy)
    if static:
        dA, dB = np.diag(T.A), np.diag(T.B)
        S.static_x = (-0.5 * dA / T.w, kT / T.w) if np.any(dA) else None
        S.static_y = (-0.5 * dB / T.w, kT / T.w) if np.any(dB) else None
    return S


def friction_scale(T: Tables, kT: float) -> Dict[str, float]:
    """K(0) of each quadratic channel and the displacement sums that set its friction."""
    _, kx = kernel(T.A, T.w, kT)
    _, ky = kernel(T.B, T.w, kT)
    w2 = T.w ** 2
    return {"K0_xQQ": float(kx.sum()), "K0_yQQ": float(ky.sum()),
            "K0_xyQ": float((T.dxy ** 2 * w2).sum()), "K0_xxQ": float((T.dxx ** 2 * w2).sum()),
            "K0_yyQ": float((T.dyy ** 2 * w2).sum())}


# ---------------------------------------------------------------------------------------------- integrator, in batches
def _ops(s, xp):
    """grad(x, y), V(x, y) and check(x, y) -> (inside, side, back) of either surface, for numpy or torch."""
    if isinstance(s, Minimal):
        from engine import coeffs, geometry, grad_poly, V_poly
        A = [float(a) for a in coeffs(s)]
        level, xc, qx, qy, nx, ny, xback = (float(v) for v in geometry(s))

        def check(x, y):
            return (x > xc) & (V_poly(A, x, y) < level), (x - qx) * nx + (y - qy) * ny, \
                (x < xback) | (abs(x) > 5.0) | (abs(y) > 5.0)
        return (lambda x, y: grad_poly(A, x, y)), (lambda x, y: V_poly(A, x, y)), check
    R, d, eps = 1.0 / s.kappa, s.d, s.eps
    level, sc = float(s.commit_level()), float(s.s_commit())
    nts2, sback = float(s.points_sn()["TS2"][1]), float(s.points_sn()["TS1"][0] - 0.4)
    atan2 = xp["atan2"]
    sqrt_ = xp["sqrt"]

    def V(x, y):
        ss, n = R * atan2(x, R - y), R - sqrt_(x * x + (R - y) ** 2)
        return 8 / 3 * ss ** 3 - 4 * ss * ss + 0.5 * n * n + (ss - d) * n * n * (n * n - 2) + eps * n

    def grad(x, y):
        rho2 = x * x + (R - y) ** 2
        rho = sqrt_(rho2)
        ss, n = R * atan2(x, R - y), R - rho
        Vs = 8 * ss * ss - 8 * ss + n * n * (n * n - 2)
        Vn = n + (ss - d) * (4 * n ** 3 - 4 * n) + eps
        return Vs * R * (R - y) / rho2 - Vn * x / rho, Vs * R * x / rho2 + Vn * (R - y) / rho

    def check(x, y):
        ss, n = R * atan2(x, R - y), R - sqrt_(x * x + (R - y) ** 2)
        V_ = 8 / 3 * ss ** 3 - 4 * ss * ss + 0.5 * n * n + (ss - d) * n * n * (n * n - 2) + eps * n
        return (ss > sc) & (V_ < level), n - nts2, (ss < sback) | (abs(n) > 3.0) | (ss > 4.0)
    return grad, V, check


class _Back:
    """numpy or torch behind the same few calls."""

    def __init__(self, backend: str, device: str = "cuda"):
        self.torch = backend == "torch"
        if self.torch:
            import torch
            self.t = torch
            self.dev = torch.device(device)
            self.fn = {"atan2": torch.atan2, "sqrt": torch.sqrt}
        else:
            self.fn = {"atan2": np.arctan2, "sqrt": np.sqrt}

    def T(self, a):
        if self.torch:
            return self.t.as_tensor(np.asarray(a, float), dtype=self.t.float64, device=self.dev)
        return np.array(a, float)

    def np(self, a):
        return a.cpu().numpy() if self.torch else np.asarray(a)

    def isbad(self, Q):
        return ~self.t.isfinite(Q).all(1) | (Q.abs() > 20.0).any(1) if self.torch else \
            ~np.isfinite(Q).all(1) | (np.abs(Q) > 20.0).any(1)


def forces_fn(T: Tables, s, back: _Back):
    grad, V, check = _ops(s, back.fn)
    w2 = back.T(T.w ** 2)
    dxx, dxy, dyy, Lx, Ly = (back.T(v) for v in (T.dxx, T.dxy, T.dyy, T.Lx, T.Ly))
    A, B = back.T(T.A), back.T(T.B)
    has_d = bool(np.any(T.dxx) or np.any(T.dxy) or np.any(T.dyy))
    has_A, has_B, has_L = bool(np.any(T.A)), bool(np.any(T.B)), bool(np.any(T.Lx) or np.any(T.Ly))
    tri = T.tri
    fx0, fy0, Kx, Ky = T.fx, T.fy, T.Kx, T.Ky

    path = T.frame == "path"
    if path:
        R0 = 1.0 / s.kappa
        at2, sq = back.fn["atan2"], back.fn["sqrt"]

    def uv(x, y):
        """The coordinates the couplings use, with d(u, v)/d(x, y) when they are the path coordinates."""
        if not path:
            return x, y, None
        rho2 = x * x + (R0 - y) ** 2
        rho = sq(rho2)
        return R0 * at2(x, R0 - y), R0 - rho, (R0 * (R0 - y) / rho2, R0 * x / rho2, -x / rho, (R0 - y) / rho)

    def force(x, y, Q, sx=0.0, sy=0.0):
        gx, gy = grad(x, y)
        u, v, J = uv(x, y)
        Fu, Fv = fx0 + sx - Kx * u, fy0 + sy - Ky * v
        X, Y = u[:, None], v[:, None]
        if has_d:
            a = 0.5 * dxx * X * X + dxy * X * Y + 0.5 * dyy * Y * Y
            dQ = Q - a
            wd = w2 * dQ
            Fu = Fu + (wd * (dxx * X + dxy * Y)).sum(1)
            Fv = Fv + (wd * (dxy * X + dyy * Y)).sum(1)
            FQ = -wd
        else:
            FQ = -w2 * Q
        if has_A:
            AQ = Q @ A
            Fu = Fu - 0.5 * (AQ * Q).sum(1)
            FQ = FQ - X * AQ
        if has_B:
            BQ = Q @ B
            Fv = Fv - 0.5 * (BQ * Q).sum(1)
            FQ = FQ - Y * BQ
        if has_L:
            Fu = Fu + (Lx * Q).sum(1)
            Fv = Fv + (Ly * Q).sum(1)
            FQ = FQ + X * Lx + Y * Ly
        for i, j, k, phi in tri:
            FQ[:, i] = FQ[:, i] - phi * Q[:, j] * Q[:, k]
            FQ[:, j] = FQ[:, j] - phi * Q[:, i] * Q[:, k]
            FQ[:, k] = FQ[:, k] - phi * Q[:, i] * Q[:, j]
        if J is None:
            return -gx + Fu, -gy + Fv, FQ
        return -gx + Fu * J[0] + Fv * J[2], -gy + Fu * J[1] + Fv * J[3], FQ

    def coupling(x, y, Q, sx=0.0, sy=0.0):
        """Everything in the potential apart from V0 and the bare harmonic energy of the modes."""
        u, v, _ = uv(x, y)
        X, Y = u[:, None], v[:, None]
        a = 0.5 * dxx * X * X + dxy * X * Y + 0.5 * dyy * Y * Y
        U = (0.5 * w2 * (a * a - 2 * a * Q)).sum(1)
        if has_A:
            U = U + 0.5 * u * ((Q @ A) * Q).sum(1)
        if has_B:
            U = U + 0.5 * v * ((Q @ B) * Q).sum(1)
        if has_L:
            U = U - u * (Lx * Q).sum(1) - v * (Ly * Q).sum(1)
        for i, j, k, phi in tri:
            U = U + phi * Q[:, i] * Q[:, j] * Q[:, k]
        return U - (fx0 + sx) * u - (fy0 + sy) * v + 0.5 * Kx * u * u + 0.5 * Ky * v * v

    force.uv = uv

    return force, coupling, V, check, w2


def sample(s, T: Tables, n: int, rng, kT: float, speed=None, occupation: str = "classical", sweeps: int = 60,
           leap: int = 20, dt: float = 0.02, backend: str = "numpy", device: str = "cuda"):
    """Start on the dividing line through TS1 with (position on the line, modes) from the Boltzmann weight of the full
    potential, by Hamiltonian Monte Carlo on that line. occupation = 'bose' or 'wigner' draws the modes independently at the
    quantum populations of each frequency (with and without the zero-point energy) instead."""
    q0, u, e, om, _ = ts_frame(s)
    M = T.M
    qw = rng.normal(0.0, sqrt(kT) / om, n)
    pu = np.full(n, float(speed)) if speed is not None else rng.rayleigh(sqrt(kT), n)
    pw = rng.normal(0.0, sqrt(kT), n)
    PQ = rng.normal(0.0, 1.0, (n, M))
    if occupation == "classical":
        Tw = np.full(M, kT)
    elif occupation == "bose":
        Tw = HBAR * T.w / np.expm1(HBAR * T.w / kT)
    elif occupation == "wigner":
        Tw = 0.5 * HBAR * T.w / np.tanh(0.5 * HBAR * T.w / kT)
    else:
        raise ValueError(occupation)
    Q = rng.normal(0.0, 1.0, (n, M)) * np.sqrt(Tw) / T.w
    PQ = PQ * np.sqrt(Tw)
    sx = sy = 0.0
    if T.static_x is not None:
        sx = (T.static_x[0] * rng.exponential(1.0, (n, T.static_x[0].size)) * T.static_x[1]).sum(1) - \
            (T.static_x[0] * T.static_x[1]).sum()
    if T.static_y is not None:
        sy = (T.static_y[0] * rng.exponential(1.0, (n, T.static_y[0].size)) * T.static_y[1]).sum(1) - \
            (T.static_y[0] * T.static_y[1]).sum()
    if occupation == "classical":
        back = _Back(backend, device)
        force, coupling, V, _, w2 = forces_fn(T, s, back)
        ex, ey = float(e[0]), float(e[1])
        q0x, q0y = float(q0[0]), float(q0[1])
        sxT, syT = (back.T(np.broadcast_to(v, (n,))) for v in (sx, sy))

        def H(q, p, Qv, P):
            x, y = q0x + q * ex, q0y + q * ey
            return 0.5 * p * p + (0.5 * P * P).sum(1) + V(x, y) + coupling(x, y, Qv, sxT, syT) + \
                (0.5 * w2 * Qv * Qv).sum(1)

        qc, Qc = back.T(qw), back.T(Q)
        acc = 0.0
        for sw in range(sweeps):
            p, P = back.T(rng.normal(0, sqrt(kT), n)), back.T(rng.normal(0, sqrt(kT), (n, M)))
            E0 = H(qc, p, Qc, P)
            q1, Q1 = qc, Qc
            h = dt * (0.8 + 0.4 * rng.uniform())

            def F(q_, Q_):
                x, y = q0x + q_ * ex, q0y + q_ * ey
                Fx, Fy, FQ = force(x, y, Q_, sxT, syT)
                return Fx * ex + Fy * ey, FQ

            f, fQ = F(q1, Q1)
            for _ in range(leap):
                p = p + 0.5 * h * f
                P = P + 0.5 * h * fQ
                q1 = q1 + h * p
                Q1 = Q1 + h * P
                f, fQ = F(q1, Q1)
                p = p + 0.5 * h * f
                P = P + 0.5 * h * fQ
            E1 = H(q1, p, Q1, P)
            dE = E1 - E0
            r = back.T(np.log(rng.uniform(size=n)))
            ok = r < -dE / kT
            if back.torch:
                ok = ok & back.t.isfinite(dE)
                qc = back.t.where(ok, q1, qc)
                Qc = back.t.where(ok[:, None], Q1, Qc)
                acc += float(ok.double().mean())
            else:
                ok = ok & np.isfinite(dE)
                qc = np.where(ok, q1, qc)
                Qc = np.where(ok[:, None], Q1, Qc)
                acc += float(ok.mean())
        qw, Q = back.np(qc), back.np(Qc)
        accept = acc / sweeps
    else:
        accept = float("nan")
    x = q0[0] + qw * e[0]
    y = q0[1] + qw * e[1]
    px = pu * u[0] + pw * e[0]
    py = pu * u[1] + pw * e[1]
    return dict(x=x, y=y, px=px, py=py, Q=Q, P=PQ, sx=sx, sy=sy, accept=accept)


def run(s, T: Tables, init: Dict, dt: float = 0.005, t_max: float = 60.0, check_every: int = 4,
        backend: str = "numpy", device: str = "cuda", frozen: bool = False) -> Dict:
    """Velocity Verlet; first-entry assignment as in ptsb.run (+1 / -1 product side, 0 returned behind TS1 or escaped,
    2 undecided, 3 a mode left the region of the surface). The energy of (x, y), of the coupling and of the modes is taken
    from the start to the moment of assignment, with the integral of p^2 + p_y^2, to measure friction directly.
    frozen = True holds every mode at its starting value: (x, y) then move on V0 + sum_j (1/2) w_j^2 (Q_j - a_j)^2, a surface
    that is not relaxed, which is the opposite limit to `run_following`."""
    back = _Back(backend, device)
    force, coupling, V, check, w2 = forces_fn(T, s, back)
    x, y, px, py, Q, P = (back.T(init[k]) for k in ("x", "y", "px", "py", "Q", "P"))
    n = x.shape[0]
    sx, sy = (back.T(np.broadcast_to(init[k], (n,))) for k in ("sx", "sy"))
    P = P * 0.0 if frozen else P
    idx = back.t.arange(n, device=back.dev) if back.torch else np.arange(n)
    out = np.full(n, 2, np.int64)
    t_out = np.full(n, np.nan)

    def Es(x, y, px, py):
        return 0.5 * (px * px + py * py) + V(x, y)

    def Em(Q, P):
        return 0.5 * (P * P + w2 * Q * Q)

    dxx, dxy, dyy = (back.T(v) for v in (T.dxx, T.dxy, T.dyy))

    def Eleft(x, y, px, py, Q, P):
        """Energy per mode in the frame of its displaced position, and the energy of following."""
        u, v, J = force.uv(x, y)
        du, dv = (px, py) if J is None else (J[0] * px + J[1] * py, J[2] * px + J[3] * py)
        X, Y, VX, VY = u[:, None], v[:, None], du[:, None], dv[:, None]
        a = 0.5 * dxx * X * X + dxy * X * Y + 0.5 * dyy * Y * Y
        ad = (dxx * X + dxy * Y) * VX + (dxy * X + dyy * Y) * VY
        return 0.5 * (P - ad) ** 2 + 0.5 * w2 * (Q - a) ** 2, 0.5 * ad * ad

    E0s = back.np(Es(x, y, px, py)).copy()
    E0c = back.np(coupling(x, y, Q, sx, sy)).copy()
    E0m = back.np(Em(Q, P)).copy()
    E0l = back.np(Eleft(x, y, px, py, Q, P)[0]).copy()
    Efs, Efc = np.full(n, np.nan), np.full(n, np.nan)
    Efm, Efl = np.full((n, T.M), np.nan), np.full((n, T.M), np.nan)
    Eff = np.full(n, np.nan)
    Ip = np.zeros(n)
    ip = back.T(np.zeros(n))
    fx, fy, fQ = force(x, y, Q, sx, sy)
    fQ = fQ * 0.0 if frozen else fQ
    for k in range(int(round(t_max / dt))):
        px = px + 0.5 * dt * fx
        py = py + 0.5 * dt * fy
        P = P + 0.5 * dt * fQ
        x = x + dt * px
        y = y + dt * py
        Q = Q + dt * P
        fx, fy, fQ = force(x, y, Q, sx, sy)
        fQ = fQ * 0.0 if frozen else fQ
        px = px + 0.5 * dt * fx
        py = py + 0.5 * dt * fy
        P = P + 0.5 * dt * fQ
        ip = ip + dt * (px * px + py * py)
        if (k + 1) % check_every == 0:
            ins, side, bk = check(x, y)
            bad = back.isbad(Q)
            stop = ins | bk | bad
            if bool(stop.any()):
                if back.torch:
                    code = back.t.where(bk, back.t.zeros_like(side), back.t.sign(side))
                    code = back.t.where(bad, 3.0 * back.t.ones_like(side), code)
                    ii = back.np(idx[stop])
                    out[ii] = back.np(code[stop]).astype(np.int64)
                else:
                    code = np.where(bad, 3.0, np.where(bk, 0.0, np.sign(side)))
                    ii = idx[stop]
                    out[ii] = code[stop].astype(np.int64)
                t_out[ii] = (k + 1) * dt
                Efs[ii] = back.np(Es(x, y, px, py)[stop])
                Efc[ii] = back.np(coupling(x, y, Q, sx, sy)[stop])
                Efm[ii] = back.np(Em(Q, P)[stop])
                el, ef = Eleft(x, y, px, py, Q, P)
                Efl[ii] = back.np(el[stop])
                Eff[ii] = back.np(ef.sum(1)[stop])
                Ip[ii] = back.np(ip[stop])
                keep = ~stop
                if not bool(keep.any()):
                    idx = idx[keep]
                    break
                x, y, px, py, Q, P, fx, fy, fQ, idx, ip, sx, sy = (
                    v[keep] for v in (x, y, px, py, Q, P, fx, fy, fQ, idx, ip, sx, sy))
    done = (out == 1) | (out == -1)
    r = {"outcome": out, "time": t_out}
    if done.any():
        # gamma_eff: energy the surface coordinates lost, against the Langevin balance d<E>/dt = -gamma (<p^2> - 2 kT)
        r["dE_surface"] = float(np.mean(Efs[done] - E0s[done]))
        r["dE_coupling"] = float(np.mean(Efc[done] - E0c[done]))
        r["dE_modes"] = float(np.mean(Efm[done].sum(1) - E0m[done].sum(1)))
        r["dE_left"] = float(np.mean(Efl[done].sum(1) - E0l[done].sum(1)))
        r["dE_left_by_mode"] = np.mean(Efl[done] - E0l[done], axis=0).tolist()
        r["E_following"] = float(np.mean(Eff[done]))
        r["int_p2"] = float(np.mean(Ip[done]))
        r["t_mean"] = float(np.mean(t_out[done]))
        Etot = (Efs[done] + Efc[done] + Efm[done].sum(1)) - (E0s[done] + E0c[done] + E0m[done].sum(1))
        r["energy_drift_rms"] = float(np.sqrt(np.mean(Etot ** 2)))
    return r


def gamma_from_energy(r: Dict, kT: float) -> float:
    """Friction read from the energy lost by the surface coordinates up to the moment of assignment, against the Langevin
    balance d<E>/dt = -gamma (<p^2> - 2 kT). It includes the reversible energy of following."""
    den = r["int_p2"] - 2 * kT * r["t_mean"]
    return float(-r["dE_surface"] / den) if den > 0 else float("nan")


def run_following(s, T: Tables, init, dt: float = 0.005, t_max: float = 60.0, check_every: int = 4,
                  backend: str = "numpy", device: str = "cuda") -> Dict:
    """The limit of modes that follow exactly: G(q) q'' = -grad V - sum_j (d_i a_j)(q'^T H_j q'), G = 1 + sum_j d a_j d a_j^T,
    H_j the constant second derivatives of a_j. No noise from the modes and no energy left in them. Classical RK4."""
    back = _Back(backend, device)
    grad, V, check = _ops(s, back.fn)
    dxx, dxy, dyy = (back.T(v) for v in (T.dxx, T.dxy, T.dyy))
    x, y, vx, vy = (back.T(v) for v in init)
    n = x.shape[0]
    idx = back.t.arange(n, device=back.dev) if back.torch else np.arange(n)
    out = np.full(n, 2, np.int64)
    t_out = np.full(n, np.nan)

    def acc(x, y, vx, vy):
        X, Y = x[:, None], y[:, None]
        ax, ay = dxx * X + dxy * Y, dxy * X + dyy * Y
        gxx, gxy, gyy = 1.0 + (ax * ax).sum(1), (ax * ay).sum(1), 1.0 + (ay * ay).sum(1)
        VX, VY = vx[:, None], vy[:, None]
        hv = dxx * VX * VX + 2.0 * dxy * VX * VY + dyy * VY * VY
        gx, gy = grad(x, y)
        fx, fy = -gx - (ax * hv).sum(1), -gy - (ay * hv).sum(1)
        det = gxx * gyy - gxy * gxy
        return (gyy * fx - gxy * fy) / det, (gxx * fy - gxy * fx) / det

    def energy(x, y, vx, vy):
        X, Y = x[:, None], y[:, None]
        ad = (dxx * X + dxy * Y) * vx[:, None] + (dxy * X + dyy * Y) * vy[:, None]
        return 0.5 * (vx * vx + vy * vy) + 0.5 * (ad * ad).sum(1) + V(x, y)

    E0 = back.np(energy(x, y, vx, vy)).copy()
    Ef = np.full(n, np.nan)
    for k in range(int(round(t_max / dt))):
        a1x, a1y = acc(x, y, vx, vy)
        a2x, a2y = acc(x + 0.5 * dt * vx, y + 0.5 * dt * vy, vx + 0.5 * dt * a1x, vy + 0.5 * dt * a1y)
        v2x, v2y = vx + 0.5 * dt * a1x, vy + 0.5 * dt * a1y
        a3x, a3y = acc(x + 0.5 * dt * v2x, y + 0.5 * dt * v2y, vx + 0.5 * dt * a2x, vy + 0.5 * dt * a2y)
        v3x, v3y = vx + 0.5 * dt * a2x, vy + 0.5 * dt * a2y
        a4x, a4y = acc(x + dt * v3x, y + dt * v3y, vx + dt * a3x, vy + dt * a3y)
        v4x, v4y = vx + dt * a3x, vy + dt * a3y
        x = x + dt / 6 * (vx + 2 * v2x + 2 * v3x + v4x)
        y = y + dt / 6 * (vy + 2 * v2y + 2 * v3y + v4y)
        vx = vx + dt / 6 * (a1x + 2 * a2x + 2 * a3x + a4x)
        vy = vy + dt / 6 * (a1y + 2 * a2y + 2 * a3y + a4y)
        if (k + 1) % check_every == 0:
            ins, side, bk = check(x, y)
            stop = ins | bk
            if bool(stop.any()):
                code = back.t.where(bk, back.t.zeros_like(side), back.t.sign(side)) if back.torch else \
                    np.where(bk, 0.0, np.sign(side))
                ii = back.np(idx[stop]) if back.torch else idx[stop]
                out[ii] = back.np(code[stop]).astype(np.int64)
                t_out[ii] = (k + 1) * dt
                Ef[ii] = back.np(energy(x, y, vx, vy)[stop])
                keep = ~stop
                if not bool(keep.any()):
                    idx = idx[keep]
                    break
                x, y, vx, vy, idx = (v[keep] for v in (x, y, vx, vy, idx))
    done = np.isfinite(Ef)
    return {"outcome": out, "time": t_out, "energy_drift_rms": float(np.sqrt(np.mean((Ef[done] - E0[done]) ** 2)))}


# ------------------------------------------------------------------------------------------------- the linear-stage law
def law(s: Minimal, T: Tables, speed: float, kT: float, t_end: float = 14.0, n_eval: int = 2801, draws: int = 4000,
        bins: int = 16, rng=None) -> Dict[str, float]:
    """The law for the cubic system, as described in the module docstring."""
    rng = rng if rng is not None else np.random.default_rng(1)
    w = T.w
    oy, ky = binned(*kernel(T.B, w, kT), bins)
    ox, kx = binned(*kernel(T.A, w, kT), bins)
    cy, cx = oy * np.sqrt(ky), ox * np.sqrt(kx)
    K0y, K0x = float(ky.sum()), float(kx.sum())
    fx0, fy0 = means_and_static(T, kT)
    jy = np.where(T.dxy != 0)[0]                   # modes of C x y Q: linear in (y, Q) along the path
    Cy = T.dxy[jy] * w[jy] ** 2
    jr = np.where((T.dxx != 0) | (T.dyy != 0))[0]  # modes displaced along the reference path
    Nr, Ny, Nb, Nx = jr.size, jy.size, oy.size, ox.size
    D = 2 + 2 * (Ny + Nb)
    wy2, wr2 = w[jy] ** 2, w[jr] ** 2
    dxx_r, dyy_r = T.dxx[jr], T.dyy[jr]
    k0 = s.hessian(0.0, 0.0)[1, 1] - K0y
    if k0 <= 0:
        return {"P_law": float("nan")}
    y0 = (fy0 - s.grad(0.0, 0.0)[1]) / k0
    vy0 = kT / k0
    g = np.r_[1.0, 0.0, np.zeros(Ny), np.zeros(Ny), (cy / oy ** 2 if Nb else np.zeros(0)), np.zeros(Nb)]
    S0 = vy0 * np.outer(g, g)
    S0[1, 1] += kT
    S0[2:2 + Ny, 2:2 + Ny] += np.diag(kT / wy2)
    S0[2 + Ny:2 + 2 * Ny, 2 + Ny:2 + 2 * Ny] += np.eye(Ny) * kT
    b0 = 2 + 2 * Ny
    S0[b0:b0 + Nb, b0:b0 + Nb] += np.diag(kT / oy ** 2) if Nb else 0
    S0[b0 + Nb:, b0 + Nb:] += np.eye(Nb) * kT
    m0 = y0 * g
    f0 = (1.0 / k0) * g                            # response to a unit static force: its start, and it is driven below

    A = np.zeros((D, D))
    A[0, 1] = 1.0
    if Ny:
        A[2:2 + Ny, 2 + Ny:2 + 2 * Ny] = np.eye(Ny)
        A[2 + Ny:2 + 2 * Ny, 2:2 + Ny] = -np.diag(wy2)
    if Nb:
        A[1, b0:b0 + Nb] = cy
        A[b0:b0 + Nb, b0 + Nb:] = np.eye(Nb)
        A[b0 + Nb:, 0] = cy
        A[b0 + Nb:, b0:b0 + Nb] = -np.diag(oy ** 2)

    o = 2 + 2 * Nr + 2 * Nx

    def rhs(t, z):
        xr, v = z[0], z[1]
        Qr, Pr = z[2:2 + Nr], z[2 + Nr:2 + 2 * Nr]
        Rx, PRx = z[2 + 2 * Nr:2 + 2 * Nr + Nx], z[2 + 2 * Nr + Nx:o]
        m, mf = z[o:o + D], z[o + D:o + 2 * D]
        S = z[o + 2 * D:].reshape(D, D)
        ar = 0.5 * dxx_r * xr * xr
        dq = Qr - ar
        Vx, Vy = s.grad(xr, 0.0)
        acc = -Vx + fx0 + K0x * xr + (wr2 * dq * dxx_r * xr).sum() + (cx * Rx).sum()
        k = s.hessian(xr, 0.0)[1, 1] - K0y - (wr2 * dq * dyy_r).sum()
        Am = A.copy()
        Am[1, 0] = -k
        if Ny:
            Cx = Cy * xr
            Am[1, 0] -= (Cx * Cx / wy2).sum()
            Am[1, 2:2 + Ny] = Cx
            Am[2 + Ny:2 + 2 * Ny, 0] = Cx
        b = np.zeros(D)
        b[1] = fy0 - Vy
        bf = np.zeros(D)
        bf[1] = 1.0
        return np.r_[v, acc, Pr, -wr2 * dq, PRx, -ox ** 2 * Rx + cx * xr, Am @ m + b, Am @ mf + bf,
                     (Am @ S + S @ Am.T).ravel()]

    def behind(t, z):
        return z[0] + 0.4
    behind.terminal = True

    z0 = np.r_[0.0, speed, np.zeros(2 * Nr), np.zeros(2 * Nx), m0, f0, S0.ravel()]
    t = np.linspace(0, t_end, n_eval)
    sol = solve_ivp(rhs, (0, t_end), z0, t_eval=t, rtol=1e-6, atol=1e-9, method="DOP853", events=(behind,))
    xr = sol.y[0]
    Qr = sol.y[2:2 + Nr]
    k_t = np.array([s.hessian(xr[i], 0.0)[1, 1] for i in range(xr.size)]) - K0y
    if Nr:
        k_t = k_t - ((Qr - 0.5 * dxx_r[:, None] * xr ** 2) * (wr2 * dyy_r)[:, None]).sum(0)
    if Ny:
        k_t = k_t + (Cy ** 2 / wy2).sum() * xr ** 2
    beyond = k_t < 0
    if not beyond.any():
        return {"P_law": float("nan")}
    k0_ = int(np.argmax(beyond))
    backs = np.where(~beyond[k0_:])[0]
    k1 = k0_ + (int(backs[0]) if backs.size else len(xr) - k0_)
    vy = sol.y[o + 2 * D]
    kb = k0_ + int(np.argmax(vy[k0_:k1]))
    my, rf = sol.y[o, kb], sol.y[o + D, kb]
    sd = sqrt(vy[kb])
    P0 = float(Phi(my / sd))
    if np.any(np.diag(T.B)):                       # the static force of each trajectory, set by the actions of the modes
        dB = np.diag(T.B)
        fs = (-0.5 * dB / w) * (rng.exponential(1.0, (draws, w.size)) * kT / w)
        fs = fs.sum(1) - fy0
        P = float(np.mean(Phi((my + fs * rf) / sd)))
    else:
        P = P0
    return {"P_law": P, "P_law_mean_force": P0, "z": float(my / sd), "x_decision": float(xr[kb]),
            "growth": float(sqrt(vy[kb] / vy0)), "k0": float(k0), "y0": float(y0), "speed_out": float(sol.y[1, kb])}


def law_thermal(s: Minimal, T: Tables, kT: float, n_speeds: int = 16, **kw) -> Dict[str, float]:
    u = (np.arange(n_speeds) + 0.5) / n_speeds
    v = np.sqrt(-2.0 * kT * np.log(1.0 - u))
    L = [law(s, T, float(vi), kT, **kw) for vi in v]
    return {"P_law": float(np.nanmean([l["P_law"] for l in L])),
            "P_law_mean_force": float(np.nanmean([l.get("P_law_mean_force", np.nan) for l in L])),
            "x_decision": float(np.nanmean([l.get("x_decision", np.nan) for l in L]))}
