"""A curved reaction valley, optionally with a pseudo-minimum before the valley-ridge inflection.

The surface is the minimal surface of ptsb.py written in curvilinear coordinates (s along the valley, n across it) and
laid along a circular arc of curvature kappa = 1/R that bends toward +n:

    X = (R - n) sin(s/R),  Y = R - (R - n) cos(s/R),     i.e. n = R - |r - (0, R)|,  s = R atan2(X, R - Y),

    V0(s, n) = 8/3 s^3 - 4 s^2 + n^2/2 + (s - d) n^2 (n^2 - 2) + eps n - a exp(-(s - s_p)^2 / (2 w^2)).

The potential depends on (s, n) exactly as the straight surface does, so the stationary points, their energies, the
intrinsic reaction coordinate and every static quantity along the path are the same as without curvature: on the line
n = 0 the force across the path is eps, and for eps = 0 the IRC runs along n = 0 to the ridge saddle TS2 and the two
product minima have the same energy. Only the kinetic energy, (1/2)[(1 - kappa n)^2 s'^2 + n'^2], carries the
curvature. The Gaussian term puts a dip in the profile along the path at s_p, before the VRI at s = 1/4 + d: a shoulder
for small a, a pseudo-minimum (an intermediate INT with an exit saddle TS1b) for larger a. A second Gaussian
+b exp(-(s - s_b)^2 / (2 w_b^2)) (with a < 0 raising the well) makes an intermediate whose exit saddle TS1b lies at the
energy of TS1, so that trajectories are held in it rather than crossing it with the energy of the descent.

Outcome +1 is the product inside the bend (n > 0, toward the centre of curvature), -1 the product outside it.

The law along the curved path. Linearizing the equations of motion about the reference s_ref(t) on n = 0, with the
fluctuations (ds, dv) along the path and (n, p) across it,

    n'' = - gamma n' - [V0_nn(s_ref, 0) - kappa^2 v^2] n - V0_n(s_ref, 0) - kappa (v^2 + C_vv) - 2 kappa v dv + noise,
    dv' = - V0_ss(s_ref, 0) ds - gamma dv - 2 kappa V0_s(s_ref, 0) n + 2 kappa v n' + noise,

where v is the reference speed and C_vv the variance of dv. The term kappa <s'^2> is the centrifugal force toward the
outside of the bend (Carpenter's dynamic matching, here as a bias of the write law); kappa^2 v^2 lowers the transverse
stiffness; 2 kappa v dv carries the fluctuations of the speed into the transverse coordinate, and the two terms of
the dv equation are the converse coupling (on the inside of the bend, n > 0, the same force along the path gives a
larger s', and a transverse velocity changes s' through the Coriolis term). In the overdamped limit
<s'^2> -> kT and the bias becomes the entropic force kappa kT of the Jacobian (1 - kappa n): the outside of the bend has
more room. The decision is read as in ptsb.linear_law.
"""
from __future__ import annotations

from math import sqrt
from typing import Dict, List, Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root

from ptsb import Phi, Surface


class CurvedValley(Surface):
    def __init__(self, kappa: float = 0.0, eps: float = 0.0, d: float = 0.25, a: float = 0.0, s_p: float = 0.25,
                 w: float = 0.06, b: float = 0.0, s_b: float = 0.38, w_b: float = 0.07):
        self.kappa, self.eps, self.d, self.a, self.s_p, self.w = (float(v) for v in (kappa, eps, d, a, s_p, w))
        self.b, self.s_b, self.w_b = float(b), float(s_b), float(w_b)
        self._points = None

    def __repr__(self):
        return f"CurvedValley(kappa={self.kappa:g}, eps={self.eps:g}, d={self.d:g}, a={self.a:g}, b={self.b:g})"

    # ---- the potential in (s, n)
    def _dip(self, s):
        return self.a * np.exp(-(s - self.s_p) ** 2 / (2 * self.w ** 2))

    def _bump(self, s):
        return self.b * np.exp(-(s - self.s_b) ** 2 / (2 * self.w_b ** 2))

    def V0(self, s, n):
        return 8 / 3 * s ** 3 - 4 * s ** 2 + 0.5 * n ** 2 + (s - self.d) * n ** 2 * (n ** 2 - 2) + self.eps * n \
            - self._dip(s) + self._bump(s)

    def V0_s(self, s, n):
        return 8 * s * s - 8 * s + n * n * (n * n - 2) + self._dip(s) * (s - self.s_p) / self.w ** 2 \
            - self._bump(s) * (s - self.s_b) / self.w_b ** 2

    def V0_n(self, s, n):
        return n + (s - self.d) * (4 * n ** 3 - 4 * n) + self.eps

    def V0_ss(self, s, n):
        return 16 * s - 8 + self._dip(s) * (1 / self.w ** 2 - (s - self.s_p) ** 2 / self.w ** 4) \
            - self._bump(s) * (1 / self.w_b ** 2 - (s - self.s_b) ** 2 / self.w_b ** 4)

    def V0_nn(self, s, n):
        return 1 + (s - self.d) * (12 * n * n - 4)

    def V0_sn(self, s, n):
        return 4 * n ** 3 - 4 * n

    # ---- the map to Cartesian coordinates
    def to_sn(self, x, y):
        if self.kappa == 0:
            return x, y
        R = 1 / self.kappa
        return R * np.arctan2(x, R - y), R - np.hypot(x, R - y)

    def to_xy(self, s, n):
        if self.kappa == 0:
            return s, n
        R = 1 / self.kappa
        return (R - n) * np.sin(s / R), R - (R - n) * np.cos(s / R)

    def V(self, x, y):
        return self.V0(*self.to_sn(x, y))

    def grad(self, x, y):
        s, n = self.to_sn(x, y)
        Vs, Vn = self.V0_s(s, n), self.V0_n(s, n)
        if self.kappa == 0:
            return Vs, Vn
        R = 1 / self.kappa
        rho2 = x * x + (R - y) ** 2
        rho = np.sqrt(rho2)
        return Vs * R * (R - y) / rho2 - Vn * x / rho, Vs * R * x / rho2 + Vn * (R - y) / rho

    def hessian(self, x, y, h: float = 1e-5):
        gxp, gyp = self.grad(x + h, y)
        gxm, gym = self.grad(x - h, y)
        hxp, hyp = self.grad(x, y + h)
        hxm, hym = self.grad(x, y - h)
        Hxx, Hxy = (gxp - gxm) / (2 * h), (hxp - hxm) / (2 * h)
        Hyx, Hyy = (gyp - gym) / (2 * h), (hyp - hym) / (2 * h)
        o = 0.5 * (Hxy + Hyx)
        return np.array([[Hxx, o], [o, Hyy]])

    # ---- stationary points, found in (s, n) and mapped to Cartesian coordinates
    def stationary_sn(self, guess):
        f = lambda q: [self.V0_s(*q), self.V0_n(*q)]
        J = lambda q: [[self.V0_ss(*q), self.V0_sn(*q)], [self.V0_sn(*q), self.V0_nn(*q)]]
        sol = root(f, np.asarray(guess, float), jac=J, tol=1e-13)
        return sol.x if sol.success and np.linalg.norm(f(sol.x)) < 1e-9 else None

    def points_sn(self) -> Dict[str, Optional[np.ndarray]]:
        if getattr(self, "_points_sn", None) is not None:
            return self._points_sn
        out = {}
        e = self.eps
        sg = np.linspace(-0.1, 0.15, 2501)
        s1 = float(sg[np.argmax(self.V0(sg, 0.0))])      # top of the profile near s = 0
        guesses = {"TS1": [s1, -e / self.V0_nn(s1, 0.0)], "TS2": [1.0, e / (4 * (1 - self.d) - 1)],
                   "P+": [1.06, 0.84], "P-": [1.06, -0.84]}
        prof = self.profile_points()
        if prof.get("INT") is not None:
            guesses["INT"] = [prof["INT"], -e / max(self.V0_nn(prof["INT"], 0), 0.1)]
            guesses["TS1b"] = [prof["TS1b"], -e / max(self.V0_nn(prof["TS1b"], 0), 0.1)]
        for k, g in guesses.items():
            q = self.stationary_sn(g)
            if q is not None:
                H = [[self.V0_ss(*q), self.V0_sn(*q)], [self.V0_sn(*q), self.V0_nn(*q)]]
                ev = np.linalg.eigvalsh(H)
                want = {"P+": "min", "P-": "min", "INT": "min"}.get(k, "saddle")
                ok = (ev > 0).all() if want == "min" else (ev[0] < 0 < ev[1])
                if k in ("P+", "P-") and np.sign(q[1]) != np.sign(g[1]):
                    ok = False
                q = q if ok else None
            out[k] = q
        self._points_sn = out
        return out

    def points(self):
        if self._points is None:
            self._points = {k: (None if q is None else np.array(self.to_xy(*q))) for k, q in self.points_sn().items()}
        return self._points

    def profile_points(self) -> Dict[str, float]:
        """Stationary points of the profile U(s) = V0(s, 0) between TS1 and the VRI: the pseudo-minimum INT and its
        exit saddle TS1b, if the dip is deep enough to make them."""
        s = np.linspace(0.02, self.vri(), 4001)
        g = self.V0_s(s, 0.0)
        idx = np.where(np.sign(g[1:]) != np.sign(g[:-1]))[0]
        out = {"INT": None, "TS1b": None}
        roots = [brentq(lambda z: self.V0_s(z, 0.0), s[i], s[i + 1]) for i in idx]
        mins = [r for r in roots if self.V0_ss(r, 0.0) > 0]
        maxs = [r for r in roots if self.V0_ss(r, 0.0) < 0 and mins and r > mins[0]]
        if mins and maxs:
            out = {"INT": mins[0], "TS1b": maxs[0]}
        return out

    def vri(self) -> float:
        return 0.25 + self.d

    # ---- product basins (in s, n)
    def commit_level(self) -> float:
        p = self.points_sn()
        vts2 = self.V0(*p["TS2"])
        return vts2 - 0.3 * (vts2 - max(self.V0(*p[k]) for k in ("P+", "P-")))

    def s_commit(self) -> float:
        return 0.5 * (self.vri() + self.points_sn()["TS2"][0])


# ---------------------------------------------------------------------------------------------- batched trajectories
def run_batch(points: List[Dict], dt: float = 0.005, check_every: int = 4, backend: str = "numpy", seed: int = 0,
              device: str = "cuda"):
    """As engine.run_batch, for curved valleys. Assignment at the first entry into a product basin (s beyond
    s_commit and V below the commit level), the side read from n against n(TS2); return behind TS1 (s < s_TS1 - 0.4)
    or leaving |n| < 3, s < 4 counts as reactant. Also returns the time of the decision."""
    sizes = [len(p["init"][0]) for p in points]
    rep = lambda f: np.concatenate([np.full(n, float(f(p))) for n, p in zip(sizes, points)])
    sf = lambda p: p["surface"]
    kap = rep(lambda p: sf(p).kappa)
    R = np.where(kap > 0, 1 / np.where(kap > 0, kap, 1.0), 1.0)
    par = {k: rep(lambda p, k=k: getattr(sf(p), k)) for k in ("eps", "d", "a", "s_p", "w", "b", "s_b", "w_b")}
    geo = {"level": rep(lambda p: sf(p).commit_level()), "sc": rep(lambda p: sf(p).s_commit()),
           "nts2": rep(lambda p: sf(p).points_sn()["TS2"][1]),
           "sback": rep(lambda p: sf(p).points_sn()["TS1"][0] - 0.4)}
    x, y, px, py = (np.concatenate([np.asarray(p["init"][i], float) for p in points]) for i in range(4))
    gamma, kT, tmax = rep(lambda p: p["gamma"]), rep(lambda p: p["kT"]), rep(lambda p: p["t_max"])
    c1 = np.exp(-gamma * dt)
    c2 = np.sqrt((1.0 - c1 * c1) * kT)
    curved = kap > 0
    steps = int(round(tmax.max() / dt))
    if backend == "torch":
        import torch
        dev = torch.device(device)
        T = lambda arr: torch.as_tensor(arr, dtype=torch.float64, device=dev)
        x, y, px, py, c1, c2, tmax, R = (T(v) for v in (x, y, px, py, c1, c2, tmax, R))
        par = {k: T(v) for k, v in par.items()}
        geo = {k: T(v) for k, v in geo.items()}
        curved = torch.as_tensor(curved, device=dev)
        gen = torch.Generator(device=dev).manual_seed(seed)
        normal = lambda: torch.randn(x.shape, generator=gen, dtype=torch.float64, device=dev)
        full = lambda v: torch.full(x.shape, v, dtype=torch.int64, device=dev)
        where, sign, exp, atan2, sqrt_ = torch.where, torch.sign, torch.exp, torch.atan2, torch.sqrt
        to_np = lambda v: v.cpu().numpy()
        as_long = lambda v: v.long()
    else:
        rng = np.random.default_rng(seed)
        normal = lambda: rng.standard_normal(x.shape)
        full = lambda v: np.full(x.shape, v, np.int64)
        where, sign, exp, atan2, sqrt_ = np.where, np.sign, np.exp, np.arctan2, np.sqrt
        to_np = lambda v: v
        as_long = lambda v: v.astype(np.int64)
    eps, d, a, sp, w, bb, sb, wb = (par[k] for k in ("eps", "d", "a", "s_p", "w", "b", "s_b", "w_b"))

    def sn(x, y):
        rho2 = x * x + (R - y) ** 2
        rho = sqrt_(rho2)
        s = where(curved, R * atan2(x, R - y), x)
        n = where(curved, R - rho, y)
        return s, n, rho, rho2

    def force(x, y):
        s, n, rho, rho2 = sn(x, y)
        dip = a * exp(-(s - sp) ** 2 / (2 * w * w))
        bump = bb * exp(-(s - sb) ** 2 / (2 * wb * wb))
        Vs = 8 * s * s - 8 * s + n * n * (n * n - 2) + dip * (s - sp) / (w * w) - bump * (s - sb) / (wb * wb)
        Vn = n + (s - d) * (4 * n ** 3 - 4 * n) + eps
        gx = where(curved, Vs * R * (R - y) / rho2 - Vn * x / rho, Vs)
        gy = where(curved, Vs * R * x / rho2 + Vn * (R - y) / rho, Vn)
        return gx, gy

    def V(s, n):
        return 8 / 3 * s ** 3 - 4 * s * s + 0.5 * n * n + (s - d) * n * n * (n * n - 2) + eps * n \
            - a * exp(-(s - sp) ** 2 / (2 * w * w)) + bb * exp(-(s - sb) ** 2 / (2 * wb * wb))

    noisy = bool((c2 > 0).any())
    out = full(2)
    tdec = full(-1)
    alive = full(1) > 0
    gx, gy = force(x, y)
    for k in range(steps):
        px = px - 0.5 * dt * gx
        py = py - 0.5 * dt * gy
        x = x + 0.5 * dt * px
        y = y + 0.5 * dt * py
        px = c1 * px
        py = c1 * py
        if noisy:
            px = px + c2 * normal() * alive
            py = py + c2 * normal() * alive
        x = x + 0.5 * dt * px
        y = y + 0.5 * dt * py
        gx, gy = force(x, y)
        px = px - 0.5 * dt * gx
        py = py - 0.5 * dt * gy
        if (k + 1) % check_every == 0:
            t = (k + 1) * dt
            s, n, _, _ = sn(x, y)
            inside = alive & (s > geo["sc"]) & (V(s, n) < geo["level"])
            sd = as_long(sign(n - geo["nts2"]))
            back = alive & ((s < geo["sback"]) | (abs(n) > 3.0) | (s > 4.0))
            out = where(back, full(0), out)
            out = where(inside, sd, out)
            tdec = where(inside | back, full(k + 1), tdec)
            expire = alive & ~back & ~inside & (t >= tmax - 1e-9)
            alive = alive & ~back & ~inside & ~expire
            if not bool(alive.any()):
                break
            px = where(alive, px, px * 0)
            py = where(alive, py, py * 0)
            gx = where(alive, gx, gx * 0)
            gy = where(alive, gy, gy * 0)
    out, tdec = to_np(out), to_np(tdec) * dt
    res, i = [], 0
    for n_ in sizes:
        res.append((out[i:i + n_], tdec[i:i + n_]))
        i += n_
    return res


# ------------------------------------------------------------------------------------- the law along a curved path
def curved_law(c: CurvedValley, speed: float, gamma: float, kT: float, start: str = "TS1", t_end: float = 12.0,
               n_eval: int = 2401, curvature: bool = True, closure: bool = True) -> Dict[str, float]:
    """The linear-stage law along the curved path (equations in the module docstring), started at the saddle `start`
    (TS1, or the exit saddle TS1b of a pseudo-minimum) with the transverse equilibrium there and the given speed.
    curvature = False drops every curvature term (the law a reaction-path analysis without the curvature coupling
    would give). closure = True adds to the mean force across the path the second-order terms -dV0_nn/ds <ds n> and
    2 kappa^2 v <dv n> in Gaussian closure: on a curved path ds and n are correlated (the inside of the bend is
    shorter, so a trajectory displaced inward is ahead along the path, where the transverse curvature is lower), and
    this correlation acts on the mean like a force toward the inside. Returns P(inside), z, and P_law = P(outside) for the reports, or NaN where the reference never passes
    the VRI (it stops in the pseudo-minimum or returns behind TS1)."""
    k = c.kappa if curvature else 0.0
    q = c.points_sn()[start]
    s0, n0 = q
    om2 = c.V0_nn(s0, n0)
    mean_n0 = -c.V0_n(s0, 0.0) / c.V0_nn(s0, 0.0)

    dVnn_ds = lambda s, h=1e-6: (c.V0_nn(s + h, 0.0) - c.V0_nn(s - h, 0.0)) / (2 * h)

    def rhs(t, z):
        s, v = z[0], z[1]
        m = z[2:6]                        # means of (ds, dv, n, p); ds and dv stay zero
        Phi_ = z[6:22].reshape(4, 4)      # propagator
        C = z[22:38].reshape(4, 4)        # covariance of the noise accumulated on the way
        Vnn = c.V0_nn(s, 0.0) - k * k * v * v
        A = np.array([[0, 1, 0, 0],
                      [-c.V0_ss(s, 0.0), -gamma, -2 * k * c.V0_s(s, 0.0), 2 * k * v],
                      [0, 0, 0, 1],
                      [0, -2 * k * v, -Vnn, -gamma]])
        Sig = C + Phi_ @ np.diag([0, var_v0, var_n0, var_p0]) @ Phi_.T
        f2 = 0.0
        if closure:   # Gaussian closure of the products ds n and dv n in the force across the path
            f2 = -dVnn_ds(s) * Sig[0, 2] + 2 * k * k * v * Sig[1, 2]
        b = np.array([0, 0, 0, -c.V0_n(s, 0.0) - k * (v * v + Sig[1, 1]) + f2])
        Q = np.diag([0, 2 * gamma * kT, 0, 2 * gamma * kT])
        return np.r_[v, -c.V0_s(s, 0.0) - gamma * v, A @ m + b, (A @ Phi_).ravel(), (A @ C + C @ A.T + Q).ravel()]

    var_v0, var_n0, var_p0 = 0.0, kT / om2, kT

    s_back = c.points_sn()["TS1"][0] - 0.4

    def behind(t, z):
        return z[0] - s_back
    behind.terminal = True

    def settled(t, z):
        return z[6 + 10] ** 2 + z[6 + 11] ** 2 - 1e24
    settled.terminal = True

    z0 = np.r_[s0, speed, 0, 0, mean_n0, 0, np.eye(4).ravel(), np.zeros(16)]
    t = np.linspace(0, t_end, n_eval)
    sol = solve_ivp(rhs, (0, t_end), z0, t_eval=t, rtol=1e-8, atol=1e-11, method="LSODA", events=(behind, settled))
    s, v = sol.y[0], sol.y[1]
    mn = sol.y[4]
    Ph = sol.y[6:22].reshape(4, 4, -1)
    C = sol.y[22:38].reshape(4, 4, -1)
    S0 = np.diag([0, var_v0, var_n0, var_p0])
    var = np.einsum("ik,kl,il->i", Ph[2].T, S0, Ph[2].T) + C[2, 2]
    amp = np.sqrt(Ph[2, 2] ** 2 + Ph[2, 3] ** 2 / om2)
    beyond = s > c.vri()
    if not beyond.any():
        return {"P_law": float("nan"), "P_in": float("nan"), "z": float("nan"), "trapped": True}
    k0 = int(np.argmax(beyond))
    back = np.where(~beyond[k0:])[0]
    k1 = k0 + (int(back[0]) if back.size else len(s) - k0)
    kb = k0 + int(np.argmax(amp[k0:k1]))
    z = mn[kb] / sqrt(var[kb])
    P_in = float(Phi(z))
    # two other readouts: where the reference first reaches the line on which trajectories are assigned
    # (s = s_commit), and where it first stops advancing (the product wall)
    kc = np.where(s[k0:k1] >= c.s_commit())[0]
    kc = k0 + int(kc[0]) if kc.size else kb
    turn = np.where((v[k0:k1 - 1] > 0) & (v[k0 + 1:k1] <= 0))[0]
    kw = k0 + int(turn[0]) + 1 if turn.size else kb
    P_commit = 1 - float(Phi(mn[kc] / sqrt(var[kc])))
    P_wall = 1 - float(Phi(mn[kw] / sqrt(var[kw])))
    return {"P_law": 1 - P_in, "P_in": P_in, "z": float(z), "mean": float(mn[kb]), "sigma": float(sqrt(var[kb])),
            "v_decision": float(v[kb]), "s_decision": float(s[kb]), "amplification": float(amp[kb]),
            "P_commit": P_commit, "P_wall": P_wall, "trapped": False}


def curved_law_thermal(c: CurvedValley, gamma: float, kT: float, n_speeds: int = 16, **kw) -> Dict[str, float]:
    """Averaged over flux-weighted thermal speeds (Rayleigh quantiles). A speed whose reference stops in the
    pseudo-minimum gives NaN; P_law is then NaN and `trapped` the share of such speeds."""
    u = (np.arange(n_speeds) + 0.5) / n_speeds
    v = np.sqrt(-2.0 * kT * np.log(1.0 - u))
    L = [curved_law(c, float(vi), gamma, kT, **kw) for vi in v]
    P = np.array([l["P_law"] for l in L])
    return {"P_law": float(P.mean()), "trapped": float(np.mean([l["trapped"] for l in L])),
            "P_commit": float(np.mean([l.get("P_commit", np.nan) for l in L])),
            "P_wall": float(np.mean([l.get("P_wall", np.nan) for l in L])),
            "v_decision": float(np.nanmean([l.get("v_decision", np.nan) for l in L]))}
