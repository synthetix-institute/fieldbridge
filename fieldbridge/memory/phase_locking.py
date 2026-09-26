"""Phase locking: the third target of co-discovery by construction.

A stable limit cycle has one flat direction, its phase (Module 6). A weak periodic modulation of the control,
p(t) = p1 + eps cos(omega_f t), acts on the phase through the phase response along the drive,
Z_p(theta) = Z(theta) . dF/dp, where Z is the gradient of the phase: the periodic solution of the adjoint equation
dZ/dt = -J^T Z, normalized by Z . F = omega_0. Averaged over the drive, only the harmonic n of Z_p with
omega_f ~ n omega_0 survives, and the phase difference psi = theta - omega_f t / n obeys the Adler equation
(Adler, Proc. IRE 34, 351 (1946))

    dpsi/dt = Delta omega - K sin phi,   phi = n psi - phi_n - pi/2,   K = eps |Z_n| / 2,
    Delta omega = omega_0 - omega_f / n,

with Z_n = a_n - i b_n the n-th Fourier coefficient of Z_p and phi_n = atan2(b_n, a_n). For |Delta omega| < K the drive
fixes the phase: psi relaxes to one of n values, 2 pi / n apart, at the rate lambda = n (K^2 - Delta omega^2)^(1/2). The
flat direction becomes a restoring one, and the phase of the drive writes the phase of the oscillator. For
|Delta omega| > K the phase slips at the frequency Omega = (Delta omega^2 - K^2)^(1/2). In canonical units
(tau = n K t, nu = Delta omega / K) the equation is dphi/dtau = nu - sin phi, and both measurements lie on one curve:
lambda / (n K) = (1 - nu^2)^(1/2) inside the locking range and Omega / K = (nu^2 - 1)^(1/2) outside.

Letters (the alphabet of codiscovery.py, acting on the same slots of the identity):

  S  a symmetry of the drift that the drive respects maps the cycle onto itself 1/m of a period later; Z_p then contains
     only the harmonics m, 2m, ..., and the drive locks at m:1 with m locked phases (m = 2: the two phases of a
     parametron)
  C  the control is moved into the range where the realization oscillates: the middle of the widest interval of
     oscillating values on a grid over the control range
  R  the state is reduced to the phase of the cycle; Z is computed from the adjoint equation, and the cycle must attract
     (a family of neutral cycles has no isolated phase)
  K  averaging over the drive gives the Adler equation; the check is the averaged phase drift measured in the driven
     realization over one period, against K cos(n psi - phi_n)
  L  the relaxation rate inside the locking range (the slowest multiplier of the stroboscopic map at the locked orbit)
     and the slip frequency outside (the time between slips of the phase difference), from simulations of the driven
     realization at two drive amplitudes; the constant is the half-width of the locking range in units of K,
     extrapolated linearly to zero amplitude

The drive amplitude is set by K = k min(omega_0, 2 kappa), k = 0.02 and 0.01, with kappa the relaxation rate of the
amplitude (the slowest transverse Floquet exponent), so that the drive is weak both against the frequency and against
the attraction of the cycle.
"""
from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from math import atan2, log, pi, sqrt
from typing import Dict, List, Optional, Tuple

import numpy as np

from . import analysis as an
from .identity import Realization

SOLVER = "LSODA"
NMAX = 8
NU_INSIDE = (-0.7, 0.0, 0.7)
NU_OUTSIDE = (-2.0, -1.3, 1.3, 2.0)
K_REL = (0.02, 0.01)
NEUTRAL = 0.99  # a transverse Floquet multiplier above this marks a family of neutral cycles
GRID = 11       # control values scanned for oscillation when the realization does not oscillate as specified


# ------------------------------------------------------------------------------------------------ integration
class _Compiled:
    """The drift of an equations realization and its exact Jacobian, compiled for one state at a time with the math
    module (single-state integration spends most of its time in the drift); network realizations use F and a
    finite-difference Jacobian."""

    def __init__(self, real: Realization):
        self.real = real
        sym = getattr(real, "symbolic", None)
        self.names = list(real.params)
        self._f = self._j = None
        if sym is not None:
            import sympy as sp
            V, P, E = sym["variables"], sym["parameters"], sym["drift"]
            self.names = [str(p) for p in P]
            self._f = sp.lambdify(V + P, E, modules="math")
            self._j = sp.lambdify(V + P, sp.Matrix(E).jacobian(V).tolist(), modules="math")

    def values(self, over: Dict) -> List[float]:
        p = {**self.real.params, **over}
        return [float(p[k]) for k in self.names]

    def f(self, q, pv) -> np.ndarray:
        if self._f is not None:
            try:
                return np.array(self._f(*q, *pv), float)
            except (ValueError, OverflowError, ZeroDivisionError, TypeError):
                pass
        return self.real.F(np.asarray(q, float), **dict(zip(self.names, pv)))

    def jac(self, q, pv) -> np.ndarray:
        if self._j is not None:
            try:
                return np.array(self._j(*q, *pv), float)
            except (ValueError, OverflowError, ZeroDivisionError, TypeError):
                pass
        return an.jacobian(self.real, np.asarray(q, float), dict(zip(self.names, pv)))


def _compiled(real: Realization) -> _Compiled:
    c = getattr(real, "_phase_compiled", None)
    if c is None:
        c = _Compiled(real)
        real._phase_compiled = c
    return c


def _solve(real: Realization, over: Dict, t_span, q0, drive=None, events=None, dense: bool = False,
           rtol: float = 1e-10, atol: Optional[float] = None):
    """solve_ivp on the realization; drive = (name, p1, eps, omega_f) modulates a parameter as p1 + eps cos(omega_f t).
    Angles are integrated without wrapping; the drift is periodic in them."""
    from scipy.integrate import solve_ivp
    atol = 1e-12 * an.base_length(real) if atol is None else atol
    c = _compiled(real)
    pv = c.values(over)
    if drive is None:
        def f(t, q):
            return c.f(q, pv)

        def jac(t, q):
            return c.jac(q, pv)
    else:
        name, p1, eps, wf = drive
        k = c.names.index(name)

        def at(t):
            v = list(pv)
            v[k] = p1 + eps * np.cos(wf * t)
            return v

        def f(t, q):
            return c.f(q, at(t))

        def jac(t, q):
            return c.jac(q, at(t))
    q0 = np.asarray(q0, float)
    if events is not None and abs(events(t_span[0], q0)) < 1e-9 * an.base_length(real):
        # a start on the section leaves the sign of the event function at t0 to rounding, and the root bracket of the
        # first step can fail: begin a short step later
        t1 = t_span[0] + 1e-6 * (t_span[1] - t_span[0])
        q0 = solve_ivp(f, (t_span[0], t1), q0, method=SOLVER, jac=jac, rtol=rtol, atol=atol).y[:, -1]
        t_span = (t1, t_span[1])
    return solve_ivp(f, t_span, q0, method=SOLVER, jac=jac, rtol=rtol, atol=atol, events=events, dense_output=dense)


def _event(j: int, level: float, torus: bool):
    """Upward crossing of the section q_j = level (for a rotating angle, of the angle level modulo 2 pi)."""
    if torus:
        def ev(t, q):
            return float(np.sin(q[j] - level))
    else:
        def ev(t, q):
            return float(q[j] - level)
    ev.direction = 1.0
    return ev


def _wrap_pi(x):
    return (np.asarray(x) + pi) % (2 * pi) - pi


# ------------------------------------------------------------------------------------------------ the cycle
def find_cycle(real: Realization, over: Dict, rng, n_starts: int = 3, chunk: float = 100.0,
               t_max: float = 3000.0) -> Optional[Dict[str, object]]:
    """A stable periodic orbit reached from random starts, or None if every start settles on a stable state or the
    motion does not become periodic. The section is the mid-level of the coordinate with the largest excursion."""
    scale = an.base_length(real)
    probe = real.carrier.sample(rng, 16)
    with np.errstate(all="ignore"):
        f_ref = float(np.nanmedian(np.linalg.norm(real.F(probe, **over), axis=-1))) or 1.0
    for _ in range(n_starts):
        q = np.asarray(real.carrier.sample(rng, 1)[0], float)
        t, sect, times, states, silent = 0.0, None, [], [], 0
        while t < t_max:
            ev = _event(*sect) if sect is not None else None
            with np.errstate(all="ignore"):
                sol = _solve(real, over, (0.0, chunk), q, events=ev, rtol=1e-9, atol=1e-11 * scale)
            if sol.status < 0 or not np.all(np.isfinite(sol.y)):
                break
            if ev is not None and sol.t_events[0].size:
                times += list(t + sol.t_events[0])
                states += [np.asarray(y, float) for y in sol.y_events[0]]
            silent = silent + 1 if ev is not None and sol.t_events[0].size == 0 else 0
            q, t = sol.y[:, -1], t + chunk
            if np.linalg.norm(real.F(q, **over)) < 1e-7 * f_ref:
                break  # settled on a stable state
            Y = sol.y[:, sol.t >= 0.5 * chunk]
            if Y.shape[1] < 8:
                continue
            span = Y.max(axis=1) - Y.min(axis=1)
            j = int(np.argmax(span))
            if span[j] < 1e-4 * scale:
                continue
            if sect is None or silent >= 2:  # no section yet, or the motion no longer crosses it
                torus = real.carrier.kind == "torus" and span[j] > real.carrier.period
                sect = (j, 0.0 if torus else float(0.5 * (Y[j].max() + Y[j].min())), torus)
                times, states, silent = [], [], 0
                continue
            if len(times) >= 6:
                P = np.diff(times[-6:])
                X = np.array(states[-6:])
                dx = max(float(real.carrier.distance(X[k], X[-1])) for k in range(len(X) - 1))
                if np.ptp(P) < 1e-6 * P[-1] and dx < 1e-5 * scale:
                    return {"period": float(P[-1]), "state": X[-1], "coordinate": sect[0], "level": sect[1],
                            "torus": bool(sect[2])}
    return None


def refine_cycle(real: Realization, over: Dict, cyc: Dict) -> Dict[str, object]:
    """Period and orbit at high accuracy; the reference point (phase 0) is the crossing of the section."""
    scale = an.base_length(real)
    ev = _event(cyc["coordinate"], cyc["level"], cyc["torus"])
    T, x0 = float(cyc["period"]), np.asarray(cyc["state"], float)
    for _ in range(2):
        sol = _solve(real, over, (0.0, 1.5 * T), x0, events=ev, rtol=1e-12, atol=1e-14 * scale)
        k = next(i for i, te in enumerate(sol.t_events[0]) if te > 0.5 * T)
        T, x0 = float(sol.t_events[0][k]), np.asarray(sol.y_events[0][k], float)
    orbit = _solve(real, over, (0.0, T), x0, dense=True, rtol=1e-12, atol=1e-14 * scale)
    return {**cyc, "period": T, "omega": 2 * pi / T, "state": x0, "orbit": orbit.sol,
            "closure": float(real.carrier.distance(orbit.y[:, -1], x0)),
            "amplitude": float(np.ptp(orbit.y[cyc["coordinate"]])) if not cyc["torus"] else float(real.carrier.period)}


def floquet(real: Realization, over: Dict, cyc: Dict) -> Dict[str, object]:
    """Multipliers of the free cycle from the monodromy matrix (central differences of the flow over one period). One
    multiplier is 1 (the phase); the largest of the others sets the relaxation rate kappa of the amplitude."""
    scale = an.base_length(real)
    x0, T = np.asarray(cyc["state"], float), float(cyc["period"])
    d, h = x0.size, 1e-5 * scale
    M = np.empty((d, d))
    for k in range(d):
        e = np.zeros(d)
        e[k] = h
        yp = _solve(real, over, (0.0, T), x0 + e, rtol=1e-12, atol=1e-14 * scale).y[:, -1]
        ym = _solve(real, over, (0.0, T), x0 - e, rtol=1e-12, atol=1e-14 * scale).y[:, -1]
        M[:, k] = (yp - ym) / (2 * h)
    mu = np.linalg.eigvals(M)
    order = np.argsort(np.abs(mu - 1.0))
    rest = np.abs(mu[order[1:]])
    other = float(rest.max()) if rest.size else 0.0
    kappa = float(-log(other) / T) if 0.0 < other < 1.0 else (float("inf") if other == 0.0 else 0.0)
    return {"multipliers": [complex(m) for m in mu[order]], "phase_multiplier": complex(mu[order[0]]),
            "largest_other": other, "kappa": kappa}


def operating_point(real: Realization, rng) -> Tuple[Optional[float], Optional[Dict], bool, List[Dict]]:
    """The control value at which the drive acts: the specified value if the realization oscillates there, otherwise
    the middle of the widest interval of oscillating values on a grid over the control range."""
    name = real.control
    p0 = float(real.params[name])
    cyc = find_cycle(real, {name: p0}, rng)
    if cyc is not None:
        return p0, cyc, False, []
    grid = np.linspace(*real.control_range, GRID)
    scan = [{"value": float(v), "oscillates": find_cycle(real, {name: float(v)}, rng) is not None} for v in grid]
    runs, cur = [], []
    for row in scan:
        if row["oscillates"]:
            cur.append(row["value"])
        elif cur:
            runs.append(cur)
            cur = []
    if cur:
        runs.append(cur)
    if not runs:
        return None, None, True, scan
    best = max(runs, key=lambda r: (len(r), -abs(np.mean(r) - p0)))
    mid = float(0.5 * (best[0] + best[-1]))
    for v in [mid] + sorted(best, key=lambda v: abs(v - mid)):
        cyc = find_cycle(real, {name: v}, rng)
        if cyc is not None:
            return float(v), cyc, True, scan
    return None, None, True, scan


# ------------------------------------------------------------------------------------------------ phase response
def phase_response(real: Realization, over: Dict, cyc: Dict, name: str, n_grid: int = 256) -> Dict[str, object]:
    """Z on the cycle from the adjoint equation integrated backward over several periods, normalized by
    Z . F = omega_0; Z_p = Z . dF/dp along the modulated parameter, and its Fourier harmonics up to NMAX."""
    from scipy.integrate import solve_ivp
    T, w0, orbit = float(cyc["period"]), float(cyc["omega"]), cyc["orbit"]
    x0 = np.asarray(cyc["state"], float)
    kap = float(cyc.get("kappa", float("inf")))
    periods = 2 if not np.isfinite(kap) else int(np.clip(np.ceil(16.0 / max(kap * T, 1e-9)), 2, 80))
    c = _compiled(real)
    pv = c.values(over)

    def jac(t, z):
        return -c.jac(orbit(t % T), pv).T

    def rhs(t, z):
        return jac(t, z) @ z

    f0 = real.F(x0, **over)
    sol = solve_ivp(rhs, (periods * T, 0.0), w0 * f0 / float(f0 @ f0), method=SOLVER, jac=jac, rtol=1e-11,
                    atol=1e-13, dense_output=True)
    ts = T * np.arange(n_grid) / n_grid
    Z = sol.sol(ts).T
    X = orbit(ts).T
    Fx = real.F(X, **over)
    norm = np.sum(Z * Fx, axis=1) / w0
    Z = Z / norm[:, None]
    p1 = float(over[name])
    dp = 1e-6 * max(1.0, abs(p1))
    dFdp = (real.F(X, **{**over, name: p1 + dp}) - real.F(X, **{**over, name: p1 - dp})) / (2 * dp)
    Zp = np.sum(Z * dFdp, axis=1)
    theta = 2 * pi * np.arange(n_grid) / n_grid
    a = np.array([2.0 * np.mean(Zp * np.cos(m * theta)) for m in range(NMAX + 1)])
    b = np.array([2.0 * np.mean(Zp * np.sin(m * theta)) for m in range(NMAX + 1)])
    return {"theta": theta, "Z": Z, "Zp": Zp, "a": a, "b": b, "abs": np.hypot(a, b),
            "normalization_spread": float(np.ptp(norm)), "backward_periods": periods}


def drive_symmetry(real: Realization, over: Dict, cyc: Dict, name: str, rng) -> Optional[Dict[str, object]]:
    """A symmetry g of the drift (permutation with signs) that holds for every value of the modulated parameter and
    maps the cycle onto itself a fraction 1/m of a period later. Then Z_p(theta + 2 pi / m) = Z_p(theta): the drive acts
    only through the harmonics m, 2m, ... The largest such m is returned."""
    from .predict import symmetries
    scale = an.base_length(real)
    p1 = float(over[name])
    found = symmetries(replace(real, params={**real.params, **over}), rng)
    p2 = p1 + 0.1 * max(1.0, abs(p1))
    qs = real.carrier.sample(rng, 12)
    T, orbit = float(cyc["period"]), cyc["orbit"]
    ts = np.linspace(0.0, T, 4001)[:-1]
    X = orbit(ts).T
    x0 = np.asarray(cyc["state"], float)
    best = None
    for g in found:
        if "permutation" not in g:
            continue
        perm, sgn = list(g["permutation"]), np.asarray(g["signs"], float)

        def act(q, perm=perm, sgn=sgn):
            gq = sgn * np.asarray(q)[..., perm]
            return real.carrier.wrap(gq) if real.carrier.kind == "torus" else gq

        F2 = real.F(qs, **{**over, name: p2})
        mismatch = np.max(np.abs(real.F(act(qs), **{**over, name: p2}) - sgn * F2[:, perm]))
        if mismatch > 1e-7 * max(1.0, np.abs(F2).max()):
            continue  # the drive breaks this symmetry
        dist = real.carrier.distance(X, act(x0)[None, :])
        k = int(np.argmin(dist))
        frac = float(ts[k] / T)
        if dist[k] > 1e-3 * scale or min(frac, 1.0 - frac) < 1e-3:
            continue  # g maps the cycle to another cycle, or onto itself without a shift in time
        m = Fraction(frac).limit_denominator(12).denominator
        if m > 1 and (best is None or m > best["m"]):
            best = {"m": int(m), "shift": frac, "permutation": perm, "signs": g["signs"], "order": g["order"]}
    return best


# ------------------------------------------------------------------------------------------------ measurements
def asymptotic_phase(real: Realization, over: Dict, cyc: Dict, q) -> float:
    """Phase of a state near the cycle: relax on the free cycle, then compare a crossing with the reference orbit."""
    T, w0 = float(cyc["period"]), float(cyc["omega"])
    n_rel = int(cyc.get("relax_periods", 1))
    ev = _event(cyc["coordinate"], cyc["level"], cyc["torus"])
    sol = _solve(real, over, (0.0, (n_rel + 1.6) * T), q, events=ev, rtol=1e-11)
    tc = sol.t_events[0]
    tc = tc[tc >= n_rel * T] if n_rel > 0 else tc
    return float((-w0 * tc[0]) % (2 * pi))


def averaged_drift(real: Realization, over: Dict, cyc: Dict, name: str, n: int, eps: float,
                   n_phases: int = 16) -> Dict[str, List[float]]:
    """The phase difference gained over one period of the oscillator (n periods of the drive at omega_f = n omega_0),
    started on the cycle at the phases psi_j: to first order in eps it is T K cos(n psi - phi_n)."""
    T, w0, orbit = float(cyc["period"]), float(cyc["omega"]), cyc["orbit"]
    p1 = float(over[name])
    psi, h = [], []
    for j in range(n_phases):
        s = (2 * pi / n) * j / n_phases
        sol = _solve(real, over, (0.0, T), orbit((s / w0) % T), drive=(name, p1, eps, n * w0), rtol=1e-11)
        gain = float(_wrap_pi(asymptotic_phase(real, over, cyc, sol.y[:, -1]) - s))
        # psi moves by the gain during the period: the average is taken at the midpoint
        psi.append(s + 0.5 * gain)
        h.append(gain / T)
    return {"psi": psi, "h": h}


def locked_rate(real: Realization, over: Dict, cyc: Dict, name: str, n: int, eps: float, K: float, phi_n: float,
                nu: float) -> Dict[str, float]:
    """Inside the locking range: the locked orbit is the fixed point of the stroboscopic map over n periods of the
    drive (Newton from the first-order locked phase); the least contracted multiplier mu gives
    lambda = -ln|mu| / (n T_drive)."""
    scale = an.base_length(real)
    T, w0, orbit = float(cyc["period"]), float(cyc["omega"]), cyc["orbit"]
    p1 = float(over[name])
    wf = n * (w0 - nu * K)
    Tm = 2 * pi * n / wf
    drive = (name, p1, eps, wf)
    psi = (phi_n + np.arccos(-nu)) / n
    q = np.asarray(orbit((psi / w0) % T), float)
    d, hq = q.size, 1e-6 * scale

    def P(y):
        return _solve(real, over, (0.0, Tm), y, drive=drive, rtol=1e-12, atol=1e-14 * scale).y[:, -1]

    def DP(y):
        cols = []
        for k in range(d):
            e = np.zeros(d)
            e[k] = hq
            cols.append((P(y + e) - P(y - e)) / (2 * hq))
        return np.array(cols).T

    res = float("inf")
    for it in range(8):
        G = real.carrier.diff(q, P(q))
        res = float(np.linalg.norm(G))
        J = DP(q)
        if res < 1e-11 * scale:
            break
        q = q + np.linalg.solve(J - np.eye(d), -G)
    mu = np.linalg.eigvals(J)
    top = float(np.max(np.abs(mu)))
    return {"nu": float(nu), "rate": float(-log(top) / Tm), "residual": res, "newton": it + 1}


def slip_frequency(real: Realization, over: Dict, cyc: Dict, name: str, n: int, eps: float, K: float,
                   nu: float) -> Dict[str, float]:
    """Outside the locking range: psi = theta - omega_f t / n is read at the crossings of the section (theta = 2 pi k)
    and interpolated; the slip frequency is 2 pi / n divided by the time between successive advances of psi by
    2 pi / n (two slips after the first)."""
    T, w0 = float(cyc["period"]), float(cyc["omega"])
    p1 = float(over[name])
    wf = n * (w0 - nu * K)
    drive = (name, p1, eps, wf)
    ev = _event(cyc["coordinate"], cyc["level"], cyc["torus"])
    step, s = 2 * pi / n, float(np.sign(nu))
    t_slip = step / (K * sqrt(nu * nu - 1.0))
    # start half a cycle away from the section (theta = pi), so that the first crossing is theta = 2 pi
    q, t, k = np.asarray(cyc["orbit"](0.5 * T), float), 0.0, 0
    psi_prev, t_prev, target, hits = pi, 0.0, pi + s * step, []
    chunk = max(t_slip / 2.0, 20 * T)
    while len(hits) < 3 and t < 8 * t_slip + 50 * T:
        sol = _solve(real, over, (t, t + chunk), q, drive=drive, events=ev, rtol=1e-10)
        for te in sol.t_events[0]:
            k += 1
            psi = 2 * pi * k - wf * te / n
            while (s > 0 and psi >= target) or (s < 0 and psi <= target):
                hits.append(t_prev + (target - psi_prev) / (psi - psi_prev) * (te - t_prev))
                target += s * step
            psi_prev, t_prev = psi, float(te)
        q, t = sol.y[:, -1], float(sol.t[-1])
    if len(hits) < 3:
        return {"nu": float(nu), "frequency": float("nan")}
    return {"nu": float(nu), "frequency": float(2 * step / (hits[2] - hits[0]))}


def adler_fit(points: List[Tuple[float, int, float]]) -> Dict[str, float]:
    """Half-width kappa and centre delta of the locking range (in units of K) from (nu, side, y): side +1 inside with
    y = lambda / (n K), (nu - delta)^2 + y^2 = kappa^2; side -1 outside with y = Omega / K, (nu - delta)^2 - y^2 =
    kappa^2. Linear in (kappa^2 - delta^2, delta) after expansion."""
    nu = np.array([p[0] for p in points])
    Y = nu ** 2 + np.array([p[1] * p[2] ** 2 for p in points])
    A = np.stack([np.ones_like(nu), 2 * nu], axis=1)
    (c0, delta), *_ = np.linalg.lstsq(A, Y, rcond=None)
    kappa = sqrt(max(c0 + delta ** 2, 0.0))
    resid = Y - A @ np.array([c0, delta])
    return {"half_width": float(kappa), "centre": float(delta), "rms": float(np.sqrt(np.mean(resid ** 2)))}


def adler_law(real: Realization, over: Dict, cyc: Dict, name: str, n: int, zn: float, phi_n: float,
              k_rel=K_REL) -> Dict[str, object]:
    """Relaxation rates inside and slip frequencies outside the locking range at two drive amplitudes; the half-width
    in units of K = eps |Z_n| / 2, extrapolated linearly in the amplitude to zero. Its uncertainty is the difference
    between the extrapolation and the value at the smaller amplitude."""
    w0, kap = float(cyc["omega"]), float(cyc.get("kappa", float("inf")))
    fits = []
    for kr in k_rel:
        K = kr * min(w0, 2.0 * kap)
        eps = 2.0 * K / zn
        pts, rows = [], []
        for nu in NU_INSIDE:
            r = locked_rate(real, over, cyc, name, n, eps, K, phi_n, nu)
            rows.append({**r, "side": "inside", "canonical": r["rate"] / (n * K)})
            pts.append((nu, 1, r["rate"] / (n * K)))
        for nu in NU_OUTSIDE:
            r = slip_frequency(real, over, cyc, name, n, eps, K, nu)
            rows.append({**r, "side": "outside", "canonical": r["frequency"] / K})
            if np.isfinite(r["frequency"]):
                pts.append((nu, -1, r["frequency"] / K))
        fits.append({"K": K, "K_over_omega": K / w0, "eps": eps, "rows": rows, **adler_fit(pts)})
    (big, small) = sorted(fits, key=lambda f: -f["eps"])[:2]
    ratio = small["eps"] / big["eps"]
    const = (small["half_width"] - ratio * big["half_width"]) / (1.0 - ratio)
    return {"fits": fits, "constant": float(const), "stderr": float(abs(const - small["half_width"])),
            "centre": float(small["centre"]), "expected": 1.0}


# ------------------------------------------------------------------------------------------------ derivation
def derive(real: Realization, rng, check_law: bool = True) -> Dict[str, object]:
    """Derive phase locking in one realization; returns the word, the steps and the certificate."""
    from .codiscovery import _field, _short, derivation_class, describe_symmetry, route_string
    word: List[str] = []
    steps: List[Dict[str, object]] = []
    notes: List[str] = []
    out: Dict[str, object] = {"name": real.name, "field": _field(real), "status": "obstructed"}

    def stop(reason: str) -> Dict[str, object]:
        out.update(obstruction=reason, obstruction_short=_short(reason))
        return {**out, "word": "".join(word), "route": route_string(steps), "class": derivation_class(steps),
                "steps": steps, "notes": notes}

    if not real.control:
        return stop("no control parameter: nothing can be modulated")
    if an.is_scale_control(real, rng):
        return stop("the control multiplies the whole drift: modulating it only changes the rate of time, so the "
                    "phase response along it is constant and has no harmonic that could fix a phase")
    name = real.control
    p1, cyc, continued, scan = operating_point(real, rng)
    if cyc is None:
        lo, hi = real.control_range
        return stop(f"no oscillation: at {GRID} values of {name} in [{lo:g}, {hi:g}] every start settles on a stable "
                    f"state or does not become periodic")
    over = {name: p1}
    if continued:
        osc = [row["value"] for row in scan if row["oscillates"]]
        word.append("C")
        steps.append({"letter": "C", "control": name, "value": float(p1),
                      "oscillating_values": [float(min(osc)), float(max(osc))]})
    cyc = refine_cycle(real, over, cyc)
    fl = floquet(real, over, cyc)
    cyc["kappa"] = fl["kappa"]
    word.append("R")
    r_step = {"letter": "R", "kind": "phase", "eliminated_directions": real.carrier.dim - 1,
              "period": cyc["period"], "omega": cyc["omega"], "largest_transverse_multiplier": fl["largest_other"],
              "kappa": fl["kappa"] if np.isfinite(fl["kappa"]) else None}
    steps.append(r_step)
    if fl["largest_other"] > NEUTRAL:
        r_step["kind"] = "neutral cycles"
        return stop("a family of neutral cycles: the transverse Floquet multiplier is 1 within the accuracy of the "
                    f"monodromy matrix (modulus {fl['largest_other']:.3f}), so the orbits next to the cycle are cycles "
                    "as well (a conserved quantity); the amplitude is as flat as the phase, and there is no isolated "
                    "phase to fix")
    cyc["relax_periods"] = 1 if not np.isfinite(fl["kappa"]) else int(np.clip(np.ceil(18.0 / max(fl["kappa"] *
                                                                                                 cyc["period"], 1e-9)),
                                                                            1, 80))
    pr = phase_response(real, over, cyc, name)
    r_step["normalization_spread"] = pr["normalization_spread"]
    harm = pr["abs"][1:]
    top = float(harm.max())
    if top <= 1e-9 * max(1.0, float(np.abs(pr["Zp"]).max())):
        return stop("the drive does not act on the phase: its phase response has no harmonic up to "
                    f"{NMAX}")
    sym = drive_symmetry(real, over, cyc, name, rng)
    m = sym["m"] if sym else 1
    n = next((k for k in range(m, NMAX + 1, m) if harm[k - 1] > 1e-4 * top), None)
    if n is None:
        return stop(f"the drive acts on the phase only through harmonics above {NMAX}")
    if sym is not None:
        word.insert(0, "S")
        meaning = describe_symmetry(real, sym)
        steps.insert(0, {"letter": "S", "meaning": meaning, "shift": sym["shift"], "ratio": m,
                         "permutation": sym["permutation"], "signs": sym["signs"], "order": sym["order"]})
        notes.append(f"{meaning} maps the cycle onto itself {Fraction(sym['shift']).limit_denominator(12)} of a period "
                     f"later and holds for every value of {name}: the drive acts through the harmonics {m}, "
                     f"{2 * m}, ...")
    zn, phi_n = float(pr["abs"][n]), float(atan2(pr["b"][n], pr["a"][n]))
    w0 = float(cyc["omega"])
    K = K_REL[-1] * min(w0, 2.0 * cyc["kappa"])
    eps = 2.0 * K / zn
    avg = averaged_drift(real, over, cyc, name, n, eps)
    s = np.asarray(avg["psi"])
    h_pred = K * np.cos(n * s - phi_n)
    dev = float(np.max(np.abs(np.asarray(avg["h"]) - h_pred)) / K)
    word.append("K")
    steps.append({"letter": "K", "ratio": n, "coupling_per_unit_drive": zn / 2.0, "phase_n": phi_n,
                  "harmonics": [float(v) for v in harm], "deviation": dev,
                  "offset": float(np.mean(avg["h"]) / K), "K": K, "eps": eps})
    stride = max(1, pr["theta"].size // 64)
    canon = {"ratio": n, "omega": w0, "period": float(cyc["period"]), "coupling_per_unit_drive": zn / 2.0,
             "phase_n": phi_n, "harmonics": [float(v) for v in harm], "deviation": dev, "K": K, "eps": eps,
             "psi": avg["psi"], "h": avg["h"],
             "phi": [float(v) for v in _wrap_pi(n * s - phi_n - pi / 2)], "g": [float(v / K) for v in avg["h"]],
             "theta": [float(v) for v in pr["theta"][::stride]], "Zp": [float(v) for v in pr["Zp"][::stride]],
             "kappa": cyc["kappa"] if np.isfinite(cyc["kappa"]) else None}
    law, const = {}, {}
    if check_law:
        law = adler_law(real, over, cyc, name, n, zn, phi_n)
        const = {"constant": law["constant"], "stderr": law["stderr"], "expected": 1.0}
        word.append("L")
        steps.append({"letter": "L", "law_constant": law["constant"], "law_stderr": law["stderr"],
                      "centre": law["centre"], "amplitudes": len(law["fits"])})
    out["status"] = "reached"
    return {**out, "word": "".join(word), "route": route_string(steps), "class": derivation_class(steps),
            "steps": steps, "canonical": canon, "law": law, "law_constant": const,
            "write_point": {"param": name, "value": float(p1)}, "notes": notes}
