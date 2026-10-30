"""Product branching after a valley-ridge inflection (VRI): model surfaces, trajectories and the linear-stage write law.

Two surfaces, masses one, a transition state TS1 at the origin, a VRI on the line y = 0 ahead of it, a ridge saddle
TS2 between two product minima:

* `Minimal(eps, c, d)`: V0 = 8/3 x^3 - 4 x^2 + y^2/2 + (x - d) y^2 (y^2 - 2), with TS1 at the origin, the VRI at
  x = 1/4 + d (V0_yy(x, 0) = 1 - 4(x - d)), TS2 at (1, 0) and, for d = 0, minima at (1.107, +-0.880); plus a tilt
  eps*y (a uniform force across the path) and a cubic c*x*y^3 (a difference in the shape of the two product valleys
  with no force on y = 0). For d < -1/4 the VRI lies behind the TS line: the origin becomes a maximum in both
  directions and TS1 splits into two saddles at y = +-[(4|d| - 1)/(4|d|)]^{1/2}, one for each product.
* `Collins(xstar)`: the model of Collins, Carpenter, Ezra and Wiggins, J. Chem. Phys. 139, 154108 (2013), eq. (2.1),
  V = c0 (x^3/3 - alpha x^2/2) + omega^2 y^2 (1 - beta x)/2 + c1 x y^4 + c2 x^2 y^2 + c3 x^2 y + c4 x y^2 + c5 x y
  + c6 x y^3 with alpha = beta = 2, c0 = 3, omega^2 = 3, and c1..c6 fixed by the lower minimum at (2.4, -1.2),
  V = -7.5 and the upper minimum at (xstar, 1.2), V = -6 (the solution reproduces their Table I).

The write law. While the displacement across the path is small it obeys, along a reference path x_ref(t) on y = 0,

    y'' = - gamma y' - V_yy(x_ref, 0) y - V_y(x_ref, 0) + sqrt(2 gamma kT) xi(t),

a linear equation: a Gaussian initial state stays Gaussian, and the side is decided by the ratio of the mean pushed by
the asymmetry (the force V_y on the line y = 0, and the offset of TS1 from that line) to the spread (initial
transverse spread at TS1 and noise accumulated on the way), P = Phi(mu / sigma) once the growing solution dominates.
This is the law P = Phi(h I1 / sqrt(sigma0^2 + 2 D I2)) of FieldBridge for a second-order (inertial) carrier. Its
overdamped normal-form limit is the Kondepudi-Nelson law; its conservative limit with the TS far before the VRI is
`airy_law`.
"""
from __future__ import annotations

from math import erf, pi, sqrt
from typing import Dict, Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root


def Phi(z):
    return 0.5 * (1.0 + np.vectorize(erf)(np.asarray(z, float) / sqrt(2.0)))


# ------------------------------------------------------------------------------------------------------ surfaces
class Surface:
    """Common machinery; subclasses define V, grad, hessian and guesses for the stationary points."""
    guesses: Dict[str, list] = {}

    def stationary(self, guess) -> Optional[np.ndarray]:
        sol = root(lambda q: np.array(self.grad(q[0], q[1])), np.asarray(guess, float),
                   jac=lambda q: self.hessian(q[0], q[1]), tol=1e-13)
        return sol.x if sol.success and np.linalg.norm(self.grad(*sol.x)) < 1e-9 else None

    def points(self) -> Dict[str, Optional[np.ndarray]]:
        if getattr(self, "_points", None) is None:
            out = {}
            for name, g in self.guesses.items():
                q = self.stationary(g)
                ev = None if q is None else np.linalg.eigvalsh(self.hessian(*q))
                want_min = name.startswith("P")
                ok = q is not None and ((ev > 0).all() if want_min else (ev[0] < 0 < ev[1]))
                if ok and want_min and np.sign(q[1]) != np.sign(g[1]):
                    ok = False
                out[name] = q if ok else None
            self._points = out
        return self._points

    def vri(self) -> float:
        """x of the valley-ridge inflection on the line y = 0 (V_yy(x, 0) = 0 between TS1 and TS2)."""
        x2 = self.points()["TS2"][0]
        return brentq(lambda x: self.hessian(x, 0.0)[1, 1], 0.0, x2)

    def ridge_normal(self):
        """TS2 and its unstable direction: the side of a point is the sign of (r - r_TS2) . n."""
        q = self.points()["TS2"]
        lam, vec = np.linalg.eigh(self.hessian(*q))
        n = vec[:, 0]
        return q, n * np.sign(n[1])

    def commit_level(self) -> float:
        """Potential below which a point beyond x_commit counts as inside a product basin: 30 % of the way from TS2
        down to the higher minimum (the basins are disconnected below the TS2 energy)."""
        p = self.points()
        vts2 = self.V(*p["TS2"])
        mins = [self.V(*p[k]) for k in ("P+", "P-") if p[k] is not None]
        return vts2 - 0.3 * (vts2 - max(mins))

    def x_commit(self) -> float:
        return 0.5 * (max(self.vri(), 0.0) + self.points()["TS2"][0])

    def vy_cubic(self, x):
        """Coefficients (highest first) of V_y(x, y) as a cubic polynomial in y."""
        ys = np.array([-1.0, -0.5, 0.5, 1.0])
        return np.polyfit(ys, np.array([self.grad(x, y)[1] for y in ys]), 3)

    def transverse_minima(self, x):
        """Local minima in y of V(x, .) at fixed x (real roots of V_y with V_yy > 0)."""
        r = np.roots(self.vy_cubic(x))
        r = np.real(r[np.abs(np.imag(r)) < 1e-9])
        return sorted(float(y) for y in r if self.hessian(x, y)[1, 1] > 0)


class Minimal(Surface):
    def __init__(self, eps: float = 0.0, c: float = 0.0, d: float = 0.0):
        self.eps, self.c, self.d = float(eps), float(c), float(d)
        if d < -0.25:
            y1 = sqrt((-4.0 * d - 1.0) / (-4.0 * d))
            lo, hi = (-y1, y1) if eps >= 0 else (y1, -y1)
            ts = {"TS1": [0.0, lo], "TS1'": [0.0, hi]}
        else:
            ts = {"TS1": [0.0, -self.eps / max(1.0 + 4.0 * d, 0.05)]}
        self.guesses = {**ts, "TS2": [1.0, self.eps / 3.0], "P+": [1.107, 0.88], "P-": [1.107, -0.88]}

    def __repr__(self):
        return f"Minimal(eps={self.eps:g}, c={self.c:g}, d={self.d:g})"

    def vri(self) -> float:
        return 0.25 + self.d

    def V(self, x, y):
        return 8.0 / 3.0 * x ** 3 - 4.0 * x ** 2 + 0.5 * y ** 2 + (x - self.d) * y ** 2 * (y ** 2 - 2.0) \
            + self.eps * y + self.c * x * y ** 3

    def grad(self, x, y):
        y2 = y * y
        return (8.0 * x * x - 8.0 * x + y2 * (y2 - 2.0) + self.c * y2 * y,
                y + (x - self.d) * (4.0 * y2 * y - 4.0 * y) + self.eps + 3.0 * self.c * x * y2)

    def hessian(self, x, y):
        Vxy = 4.0 * y ** 3 - 4.0 * y + 3.0 * self.c * y ** 2
        return np.array([[16.0 * x - 8.0, Vxy],
                         [Vxy, 1.0 + (x - self.d) * (12.0 * y ** 2 - 4.0) + 6.0 * self.c * x * y]])


class Collins(Surface):
    c0, alpha, beta, om2 = 3.0, 2.0, 2.0, 3.0

    def __init__(self, xstar: float = 2.05):
        self.xstar = float(xstar)
        A, b = [], []
        for (x, y, E) in [(2.4, -1.2, -7.5), (self.xstar, 1.2, -6.0)]:
            base = self.c0 * (x ** 3 / 3 - self.alpha * x ** 2 / 2) + self.om2 / 2 * y ** 2 * (1 - self.beta * x)
            bx = self.c0 * (x ** 2 - self.alpha * x) - self.om2 / 2 * y ** 2 * self.beta
            by = self.om2 * y * (1 - self.beta * x)
            A += [[x * y ** 4, x ** 2 * y ** 2, x ** 2 * y, x * y ** 2, x * y, x * y ** 3],
                  [y ** 4, 2 * x * y ** 2, 2 * x * y, y ** 2, y, y ** 3],
                  [4 * x * y ** 3, 2 * x ** 2 * y, x ** 2, 2 * x * y, x, 3 * x * y ** 2]]
            b += [E - base, -bx, -by]
        self.k = np.linalg.solve(np.array(A), np.array(b))
        self.guesses = {"TS1": [0.0, 0.0], "TS2": [1.98, 0.09], "P+": [self.xstar, 1.2], "P-": [2.4, -1.2]}

    def __repr__(self):
        return f"Collins(xstar={self.xstar:g})"

    def V(self, x, y):
        c1, c2, c3, c4, c5, c6 = self.k
        return self.c0 * (x ** 3 / 3 - self.alpha * x ** 2 / 2) + self.om2 * y ** 2 * (1 - self.beta * x) / 2 \
            + c1 * x * y ** 4 + c2 * x ** 2 * y ** 2 + c3 * x ** 2 * y + c4 * x * y ** 2 + c5 * x * y + c6 * x * y ** 3

    def grad(self, x, y):
        c1, c2, c3, c4, c5, c6 = self.k
        Vx = self.c0 * (x ** 2 - self.alpha * x) - self.om2 * self.beta * y ** 2 / 2 + c1 * y ** 4 \
            + 2 * c2 * x * y ** 2 + 2 * c3 * x * y + c4 * y ** 2 + c5 * y + c6 * y ** 3
        Vy = self.om2 * y * (1 - self.beta * x) + 4 * c1 * x * y ** 3 + 2 * c2 * x ** 2 * y + c3 * x ** 2 \
            + 2 * c4 * x * y + c5 * x + 3 * c6 * x * y ** 2
        return Vx, Vy

    def hessian(self, x, y):
        c1, c2, c3, c4, c5, c6 = self.k
        Vxx = self.c0 * (2 * x - self.alpha) + 2 * c2 * y ** 2 + 2 * c3 * y
        Vxy = -self.om2 * self.beta * y + 4 * c1 * y ** 3 + 4 * c2 * x * y + 2 * c3 * x + 2 * c4 * y + c5 \
            + 3 * c6 * y ** 2
        Vyy = self.om2 * (1 - self.beta * x) + 12 * c1 * x * y ** 2 + 2 * c2 * x ** 2 + 2 * c4 * x + 6 * c6 * x * y
        return np.array([[Vxx, Vxy], [Vxy, Vyy]])


# --------------------------------------------------------------------------------------------- initial conditions
def ts_frame(s: Surface):
    """TS1, its unstable direction u (toward +x), stable direction w (toward +y), frequencies omega_w, |omega_u|."""
    q = s.points()["TS1"]
    lam, vec = np.linalg.eigh(s.hessian(*q))
    u, w = vec[:, 0] * np.sign(vec[0, 0]), vec[:, 1] * np.sign(vec[1, 1])
    return q, u, w, sqrt(lam[1]), sqrt(-lam[0])


def transverse_energy(kT: float, om: float, transverse: str, hbar: float) -> float:
    """Mean energy of the transverse mode at TS1 (the variance scale of the initial state)."""
    if transverse == "thermal":
        return kT
    if transverse == "wigner":
        return 0.5 * hbar * om / np.tanh(0.5 * hbar * om / kT) if kT > 0 else 0.5 * hbar * om
    return hbar * om * (0.5 + (1.0 / np.expm1(hbar * om / kT) if kT > 0 else 0.0))


def sample_ts(s: Surface, n: int, rng, kT: float, speed=None, transverse: str = "thermal", hbar: float = 0.0):
    """Initial conditions on the dividing line through TS1 normal to its unstable direction.

    speed: fixed momentum along the unstable direction, or None for the flux-weighted thermal distribution
    (Rayleigh with scale kT^{1/2}).
    transverse:
      thermal  classical Boltzmann: q ~ N(0, kT/omega^2), p ~ N(0, kT)
      wigner   Gaussian with the quantum thermal variance (hbar/2 omega) coth(hbar omega/2kT)
      qct      quasiclassical normal-mode sampling: E = (n + 1/2) hbar omega with n Boltzmann distributed and a
               uniform phase (at kT = 0 every trajectory carries exactly the zero-point energy)
    """
    q0, u, w, om, _ = ts_frame(s)
    if transverse in ("thermal", "wigner"):
        e = transverse_energy(kT, om, transverse, hbar)
        qw, pw = rng.normal(0.0, sqrt(e) / om, n), rng.normal(0.0, sqrt(e), n)
    elif transverse == "qct":
        nq = rng.geometric(1.0 - np.exp(-hbar * om / kT), n) - 1 if kT > 0 else np.zeros(n)
        E = (nq + 0.5) * hbar * om
        ph = rng.uniform(0.0, 2.0 * pi, n)
        qw, pw = np.sqrt(2.0 * E) / om * np.cos(ph), -np.sqrt(2.0 * E) * np.sin(ph)
    else:
        raise ValueError(transverse)
    pu = np.full(n, float(speed)) if speed is not None else rng.rayleigh(sqrt(kT), n)
    return (q0[0] + qw * w[0], q0[1] + qw * w[1], pu * u[0] + pw * w[0], pu * u[1] + pw * w[1])


def sample_dividing_line(s: Surface, n: int, rng, kT: float, y_span: float = 2.0, grid: int = 4001):
    """Canonical flux ensemble on the line x = 0 through the TS region: y from exp(-V(0, y)/kT) (the full,
    anharmonic weight, so that it holds on both sides of the cusp where TS1 splits into two saddles), p_y ~ N(0, kT),
    p_x Rayleigh(kT^{1/2}) > 0."""
    y = np.linspace(-y_span, y_span, grid)
    w = np.exp(-(s.V(0.0 * y, y) - s.V(0.0 * y, y).min()) / kT)
    cdf = np.cumsum(w)
    cdf /= cdf[-1]
    ys = np.interp(rng.uniform(0, 1, n), cdf, y)
    return (np.zeros(n), ys, rng.rayleigh(sqrt(kT), n), rng.normal(0.0, sqrt(kT), n))


def tst_split(s: Surface, kT: float, y_span: float = 2.0, grid: int = 4001) -> Dict[str, float]:
    """Static (transition-state) prediction when TS1 has split into one saddle per product: the flux through each
    side of the line x = 0, divided at the maximum of V(0, y) between the two saddles. Returns P(+ side) and the
    split point, or NaN where V(0, y) has a single minimum (the bifurcation lies after the TS line)."""
    y = np.linspace(-y_span, y_span, grid)
    V = s.V(0.0 * y, y)
    mins = np.where((V[1:-1] < V[:-2]) & (V[1:-1] < V[2:]))[0] + 1
    if mins.size < 2:
        return {"P_tst": float("nan"), "y_split": float("nan")}
    lo, hi = mins[y[mins] < 0], mins[y[mins] > 0]
    if lo.size == 0 or hi.size == 0:
        return {"P_tst": float("nan"), "y_split": float("nan")}
    a, b = int(lo.max()), int(hi.min())
    k = a + int(np.argmax(V[a:b + 1]))
    w = np.exp(-(V - V.min()) / kT)
    return {"P_tst": float(w[k:].sum() / w.sum()), "y_split": float(y[k])}


# --------------------------------------------------------------------------------------------------- trajectories
def run(s: Surface, init, gamma: float = 0.0, kT: float = 0.0, dt: float = 0.005, t_max: float = 40.0, rng=None,
        rule: str = "first", check_every: int = 4) -> Dict[str, np.ndarray]:
    """BAOAB Langevin (velocity Verlet at gamma = 0; kT = 0 with gamma > 0 is pure damping, as in Collins et al.).

    rule = "first": a trajectory is assigned at its first entry into a product basin (V < commit_level beyond
    x_commit), by the side of the TS2 ridge; it is "reactant" (0) if it returns behind TS1 (x < x_TS1 - 0.4) first
    or leaves |x|, |y| < 5 (the open region behind TS1, where both surfaces fall without bound).
    rule = "final": every trajectory runs to t_max and is assigned by the basin it occupies then (2 if in neither);
    returns to the reactant side still stop the trajectory. `first` always holds the first-entry assignment.
    Outcome codes: +1 / -1 product side, 0 reactant, 2 undecided."""
    x, y, px, py = (np.array(a, float) for a in init)
    n = x.size
    level, xc = s.commit_level(), s.x_commit()
    qr, nr = s.ridge_normal()
    x_back = s.points()["TS1"][0] - 0.4
    out = np.full(n, 2, int)
    first = np.full(n, 2, int)
    t_out = np.full(n, np.nan)
    idx = np.arange(n)
    c1 = np.exp(-gamma * dt)
    c2 = np.sqrt((1.0 - c1 * c1) * kT) if gamma > 0 else 0.0
    gx, gy = s.grad(x, y)
    for k in range(int(round(t_max / dt))):
        px -= 0.5 * dt * gx
        py -= 0.5 * dt * gy
        x += 0.5 * dt * px
        y += 0.5 * dt * py
        if gamma > 0:
            px *= c1
            py *= c1
            if c2 > 0:
                px += c2 * rng.standard_normal(px.size)
                py += c2 * rng.standard_normal(py.size)
        x += 0.5 * dt * px
        y += 0.5 * dt * py
        gx, gy = s.grad(x, y)
        px -= 0.5 * dt * gx
        py -= 0.5 * dt * gy
        if (k + 1) % check_every == 0:
            inside = (x > xc) & (s.V(x, y) < level)
            sd = np.sign((x - qr[0]) * nr[0] + (y - qr[1]) * nr[1]).astype(int)
            newly = inside & (first[idx] == 2)
            if newly.any():
                first[idx[newly]] = sd[newly]
                t_out[idx[newly]] = (k + 1) * dt
            back = (x < x_back) | (np.abs(x) > 5.0) | (np.abs(y) > 5.0)
            stop = back | (inside if rule == "first" else np.zeros_like(back))
            if stop.any():
                ii = idx[stop]
                out[ii] = np.where(back[stop], 0, sd[stop])
                first[ii] = np.where(back[stop] & (first[ii] == 2), 0, first[ii])
                keep = ~stop
                x, y, px, py, gx, gy, idx = x[keep], y[keep], px[keep], py[keep], gx[keep], gy[keep], idx[keep]
                if idx.size == 0:
                    break
    if idx.size and rule == "final":
        inside = (x > xc) & (s.V(x, y) < level)
        sd = np.sign((x - qr[0]) * nr[0] + (y - qr[1]) * nr[1]).astype(int)
        out[idx] = np.where(inside, sd, 2)
    return {"outcome": out, "first": first, "time": t_out}


def ratio(outcome: np.ndarray) -> Dict[str, float]:
    n = int(outcome.size)
    plus, minus = int(np.sum(outcome == 1)), int(np.sum(outcome == -1))
    m = plus + minus
    P = plus / m if m else float("nan")
    return {"n": n, "plus": plus, "minus": minus, "reactant": int(np.sum(outcome == 0)),
            "undecided": int(np.sum(outcome == 2)), "P_plus": P, "stderr": sqrt(P * (1 - P) / m) if m else float("nan")}


# -------------------------------------------------------------------------------------- the linear-stage write law
def linear_law(s: Surface, speed: float, gamma: float, kT: float, transverse: str = "thermal", hbar: float = 0.0,
               t_end: float = 12.0, n: int = 2401, initial=None) -> Dict[str, object]:
    """Mean and spread of the transverse displacement from the linear equation along the reference path.

    The reference x_ref(t) runs on the line y = 0 from x = 0 with velocity `speed` and deterministic friction gamma.
    Integrated with it: the mean (y, p) pushed by the force V_y(x_ref, 0) from the transverse equilibrium at TS1,
    y0 = -V_y(0,0)/V_yy(0,0); the fundamental matrix M of the transverse equation (response to the initial state,
    whose variances are e/omega^2 and e with e from `transverse_energy`); and the covariance S_n of the noise
    accumulated on the way. The decision is read where the growing solution has amplified most while the reference
    is beyond the VRI (the first excursion: a conservative reference returns over the VRI, after which the linear
    picture no longer selects a side). Also returned: G, the log of the amplification reached when the reference
    first stops advancing (its first turning point, i.e. the product wall) or at the end if it never turns: the
    transverse growth achieved before the motion along the path meets the wall.

    initial = (mean_y, var_y, var_p) replaces the harmonic initial state (used with `sample_dividing_line`, whose
    anharmonic Boltzmann weight near the cusp is far from the harmonic one)."""
    _, _, _, om, _ = ts_frame(s)
    om2 = s.hessian(0.0, 0.0)[1, 1]
    e = transverse_energy(kT, om, transverse, hbar)

    def rhs(t, z):
        x, v, my, mp, m00, m01, m10, m11, s00, s01, s11 = z
        Vyy = s.hessian(x, 0.0)[1, 1]
        Vx, f = s.grad(x, 0.0)
        return [v, -Vx - gamma * v, mp, -Vyy * my - gamma * mp - f,
                m10, m11, -Vyy * m00 - gamma * m10, -Vyy * m01 - gamma * m11,
                2 * s01, s11 - Vyy * s00 - gamma * s01, -2 * Vyy * s01 - 2 * gamma * s11 + 2 * gamma * kT]

    def behind(t, z):          # the reference returned behind TS1 (the profile then falls without bound)
        return z[0] + 0.4
    behind.terminal = True

    def settled(t, z):         # amplification 10^12: the side was chosen long before
        return z[4] ** 2 + z[5] ** 2 - 1e24
    settled.terminal = True

    f0 = s.grad(0.0, 0.0)[1]
    t = np.linspace(0.0, t_end, n)
    y0, vy0, vp0 = (-f0 / om2, e / om2, e) if initial is None else initial
    sol = solve_ivp(rhs, (0.0, t_end), [0.0, speed, y0, 0.0, 1, 0, 0, 1, 0, 0, 0], t_eval=t,
                    rtol=1e-9, atol=1e-12, method="LSODA", events=(behind, settled))
    t = sol.t
    x, v, my, mp, m00, m01, m10, m11, s00, s01, s11 = sol.y
    var = m00 ** 2 * vy0 + m01 ** 2 * vp0 + s00
    amp = np.sqrt(m00 ** 2 + m01 ** 2 / om2)          # growth of the y-response to a unit-energy initial state
    beyond = x > s.vri()
    if not beyond.any():
        return {"P_law": float("nan"), "z": float("nan")}
    k0 = int(np.argmax(beyond))
    back = np.where(~beyond[k0:])[0]
    k1 = k0 + (int(back[0]) if back.size else len(x) - k0)
    kb = k0 + int(np.argmax(amp[k0:k1]))
    turn = np.where((v[k0:k1 - 1] > 0) & (v[k0 + 1:k1] <= 0))[0]
    kg = k0 + int(turn[0]) + 1 if turn.size else k1 - 1
    z = my[kb] / sqrt(var[kb])
    return {"P_law": float(Phi(z)), "z": float(z), "mean": float(my[kb]), "sigma": float(sqrt(var[kb])),
            "sigma_noise": float(sqrt(max(s00[kb], 0.0))), "t_decision": float(t[kb]), "x_decision": float(x[kb]),
            "amplification": float(amp[kb]), "G": float(np.log(amp[kg])), "a": float(m00[kb]), "b": float(m01[kb]),
            "e": e, "omega_ts": om, "trace": np.c_[t, x, my, np.sqrt(var), amp][k0:k1]}


def law_thermal_speed(s: Surface, gamma: float, kT: float, n_speeds: int = 16, **kw) -> Dict[str, float]:
    """The law averaged over the flux-weighted thermal speeds (Rayleigh, scale kT^{1/2}) at the dividing line."""
    u = (np.arange(n_speeds) + 0.5) / n_speeds
    v = np.sqrt(-2.0 * kT * np.log(1.0 - u))        # Rayleigh quantiles
    laws = [linear_law(s, float(vi), gamma, kT, **kw) for vi in v]
    return {"P_law": float(np.mean([L["P_law"] for L in laws])), "G": float(np.mean([L["G"] for L in laws]))}


def dividing_line_moments(s: Surface, kT: float, y_span: float = 2.0, grid: int = 4001):
    """Mean and variance of y under exp(-V(0, y)/kT) on the line x = 0, and the momentum variance kT."""
    y = np.linspace(-y_span, y_span, grid)
    w = np.exp(-(s.V(0.0 * y, y) - s.V(0.0 * y, y).min()) / kT)
    w /= w.sum()
    m = float((w * y).sum())
    return m, float((w * (y - m) ** 2).sum()), kT


def law_dividing_line(s: Surface, gamma: float, kT: float, n_speeds: int = 16, **kw) -> Dict[str, float]:
    """The law for `sample_dividing_line`: the anharmonic moments on x = 0, averaged over Rayleigh speeds. Defined
    only while the VRI lies after the TS line (one transverse minimum at x = 0)."""
    if s.vri() <= 0:
        return {"P_law": float("nan"), "G": float("nan")}
    init = dividing_line_moments(s, kT)
    u = (np.arange(n_speeds) + 0.5) / n_speeds
    v = np.sqrt(-2.0 * kT * np.log(1.0 - u))
    laws = [linear_law(s, float(vi), gamma, kT, initial=init, **kw) for vi in v]
    return {"P_law": float(np.mean([L["P_law"] for L in laws])), "G": float(np.mean([L["G"] for L in laws]))}


def qct_probability(law: Dict, hbar: float, kT: float) -> float:
    """P(y > 0) for quasiclassical normal-mode sampling with no noise after the start. The decision variable is
    m + a q + b p with (q, p) on the ellipse of energy E = (n + 1/2) hbar omega and a uniform phase, i.e.
    m + A cos(phi) with A = (2E)^{1/2} (a^2/omega^2 + b^2)^{1/2}; summed over the Boltzmann distribution of n.
    An arcsine law: the selectivity is complete once |m| exceeds A."""
    om, a, b, m = law["omega_ts"], law["a"], law["b"], law["mean"]
    k = sqrt(2.0 * (a * a / om ** 2 + b * b))
    ns = np.arange((int(40 * kT / (hbar * om)) + 40) if kT > 0 else 1)
    w = np.exp(-ns * hbar * om / kT) if kT > 0 else np.array([1.0])
    w = w / w.sum()
    return float(sum(wq * (0.5 + np.arcsin(np.clip(m / (k * sqrt((nq + 0.5) * hbar * om)), -1, 1)) / pi)
                     for nq, wq in zip(ns, w)))


# ------------------------------------------------------------------------------------- local (normal-form) laws
def kn_local(s: Surface, gamma: float, kT: float) -> float:
    """Kondepudi-Nelson law from the local quantities at the VRI in the overdamped limit: growth rate
    lambda = -V_yy/gamma rises at r = -dV_yy/dx * xdot / gamma with xdot = -V_x/gamma; bias h = -V_y/gamma; noise
    D = kT/gamma. P(+ side) = Phi(pi^{1/4} h / (D^{1/2} r^{1/4}))."""
    xv = s.vri()
    d = 1e-5
    dVyy = (s.hessian(xv + d, 0)[1, 1] - s.hessian(xv - d, 0)[1, 1]) / (2 * d)
    xdot = -s.grad(xv, 0.0)[0] / gamma
    r = -dVyy * xdot / gamma
    h = -s.grad(xv, 0.0)[1] / gamma
    return float(Phi(pi ** 0.25 * h / (sqrt(kT / gamma) * r ** 0.25)))


def airy_law(f: float, J: float, R: float) -> float:
    """Conservative passage through a linear sign change of the transverse stiffness, omega^2(t) = -R t, entered
    adiabatically with mean transverse action J (thermal: J = kT/omega_TS): the growing (Bi) amplitude is
    -pi f R^{-2/3} from the force f = V_y and Gaussian with variance pi J R^{-1/3} from the initial oscillation, so
    P(+ side) = Phi(-pi^{1/2} f / (J R)^{1/2})."""
    return float(Phi(-sqrt(pi) * f / sqrt(J * R)))


# --------------------------------------------------------------------------------- static structure of the surface
def string_path(s: Surface, a, b, images: int = 61, iters: int = 40000, step: float = 2e-3, tol: float = 1e-6,
                start=None):
    """Zero-temperature string method (E, Ren and Vanden-Eijnden) between fixed end points a and b: images move along
    the force and are redistributed at equal arc length. `start` is an optional initial path (images x 2)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    P = np.array(start, float) if start is not None else a + np.linspace(0, 1, images)[:, None] * (b - a)
    images = len(P)
    it = 0
    for it in range(iters):
        g = np.array(s.grad(P[1:-1, 0], P[1:-1, 1])).T
        P_new = P.copy()
        P_new[1:-1] -= step * g
        d = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P_new, axis=0), axis=1))]
        L = np.linspace(0, d[-1], images)
        P_new = np.c_[np.interp(L, d, P_new[:, 0]), np.interp(L, d, P_new[:, 1])]
        moved = np.max(np.abs(P_new - P))
        P = P_new
        if moved < tol * step:
            break
    return P, s.V(P[:, 0], P[:, 1]), it


def irc(s: Surface, sign: float = 1.0, ds: float = 1e-3, n: int = 40000):
    """Steepest-descent path from TS1 along its unstable direction (toward +x for sign = +1)."""
    q, u, *_ = ts_frame(s)
    r = q + sign * 1e-3 * u
    path = [r.copy()]
    for _ in range(n):
        g = np.array(s.grad(*r))
        ng = np.linalg.norm(g)
        if ng < 1e-7:
            break
        r = r - ds * g / ng
        path.append(r.copy())
    return np.array(path)


def folds(s: Surface, x_lo: float, x_hi: float, n: int = 4000):
    """Points along x where the number of transverse minima of V(x, .) increases: the pitchfork (VRI) of a symmetric
    surface, or the fold at which a minor valley is born on an asymmetric one. Returns (x, y_new, V) for each birth,
    y_new being the new valley (the one farther from the valleys that already existed)."""
    xs = np.linspace(x_lo, x_hi, n)
    counts = [len(s.transverse_minima(x)) for x in xs]
    out = []
    for i in range(1, n):
        if counts[i] > counts[i - 1]:
            new, old = s.transverse_minima(xs[i]), s.transverse_minima(xs[i - 1])
            y_new = max(new, key=lambda y: min(abs(y - o) for o in old)) if old else new[0]
            out.append((float(xs[i]), float(y_new), float(s.V(xs[i], y_new))))
    return out
