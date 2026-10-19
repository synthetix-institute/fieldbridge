"""Two outcomes of a material without a stored state: chaotic loss of the initial state, and blow-up in finite time.

Chaotic loss of the initial state. On a chaotic attractor two states that differ by delta separate as
delta e^{lambda_1 t}, with lambda_1 > 0 the largest Lyapunov exponent. A preparation known to the resolution delta is
therefore distinguishable from its neighbours only until the separation reaches the size of the attractor, after the
predictability time

    T(delta) = (1/lambda_1) ln(Delta / delta)                                         (Lorenz, 1969)

The knowledge of the initial state is lost at the rate lambda_1 / ln 2 bits per unit time, without noise. The test:
lambda_1 is measured from renormalized pairs (Benettin et al., 1980), T(delta) from separate pairs started on the
attractor at several delta, and the slope of T against ln(1/delta) times lambda_1 must be 1.

Blow-up in finite time. A state of a polynomial drift that escapes to infinity in finite time grows near the blow-up
time t* as x_j ~ c_j (t* - t)^(-p_j). Inserting this into dx_i/dt = sum_k a_ik prod_j x_j^(m_kj) gives the dominant
balance (Kovalevskaya, 1889; Goriely and Hyde, 2000): in every equation the largest power among the monomials equals
p_i + 1,

    max_k (m_k . p) = p_i + 1,   p_i c_i = sum over the dominant monomials of a_ik prod_j c_j^(m_kj),

with a real solution c. The growth exponent of the state is p = max_j p_j, so along the escape
d ln|F| / d ln|x| = 1 + 1/p, and the escape time from a start scaled by s falls as s^(-1/p). The total degree n of the
drift predicts p = 1/(n - 1) only when the escape runs along the terms of highest degree; x'' = x^2 written as
x' = y, y' = x^2 has degree 2 but p = 3 (x ~ (t* - t)^(-2), y ~ (t* - t)^(-3)). The test: the balances of the drift
are found by enumerating, in every equation, the monomial that dominates, solving the linear equations for p and the
nonlinear ones for c; the measured growth exponent of escaping states must equal one of them, and the escape times
of scaled starts must scale as predicted.

Limits. Polynomial drifts only for the balances (a drift with other functions is measured but not predicted). The
attractor test is numerical: a chaotic label needs lambda_1 above its spread between starts, and a growth exponent
needs the escape to reach amplitudes where the dominant terms rule.
"""
from __future__ import annotations

import itertools
from typing import Dict, List, Optional, Sequence

import numpy as np

from .identity import Realization


# ------------------------------------------------------------------------------------------------ integration
def _rk4(F, x, dt):
    k1 = F(x); k2 = F(x + 0.5 * dt * k1); k3 = F(x + 0.5 * dt * k2); k4 = F(x + dt * k3)
    return x + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)


def _F(real: Realization, over: Optional[Dict[str, float]] = None):
    over = over or {}
    return lambda x: real.F(x, **over)


# ------------------------------------------------------------------------------------------------ chaos
def settle(real: Realization, starts: np.ndarray, T: float, dt: float, over=None) -> np.ndarray:
    """Integrate the starts (rows) for a time T; rows that escape (|x| > 1e8 or not finite) are dropped."""
    F = _F(real, over)
    x = np.array(starts, float)
    for _ in range(int(T / dt)):
        x = _rk4(F, x, dt)
        bad = ~np.all(np.isfinite(x), axis=1) | (np.abs(x).max(axis=1) > 1e8)
        if bad.any():
            x = x[~bad]
            if not len(x):
                break
    return x


def lyapunov(real: Realization, starts: np.ndarray, T: float, dt: float, d0: float = 1e-8, every: int = 10,
             over=None, seed: int = 1) -> Dict[str, object]:
    """Largest Lyapunov exponent of each start (rows) by renormalized neighbours, with the value over the first half."""
    F = _F(real, over)
    rng = np.random.default_rng(seed)
    x = np.array(starts, float)
    u = rng.standard_normal(x.shape)
    u *= d0 / np.linalg.norm(u, axis=1, keepdims=True)
    y = x + u
    s, half, n = np.zeros(len(x)), None, int(T / dt)
    for k in range(1, n + 1):
        x, y = _rk4(F, x, dt), _rk4(F, y, dt)
        if k % every == 0:
            d = np.linalg.norm(y - x, axis=1)
            if not np.all(np.isfinite(d)) or np.any(d == 0):
                return {"lambda1": None}
            s += np.log(d / d0)
            y = x + (y - x) * (d0 / d)[:, None]
            if k == (n // (2 * every)) * every:
                half = s / (k * dt)
    lam = s / (n * dt)
    return {"lambda1": lam.tolist(), "first_half": None if half is None else half.tolist(),
            "median": float(np.median(lam)), "spread": float(lam.max() - lam.min())}


def predictability(real: Realization, points: np.ndarray, lam: float, size: float, dt: float,
                   deltas: Sequence[float] = (1e-9, 1e-7, 1e-5, 1e-3), reach: float = 0.1, over=None,
                   seed: int = 2) -> Dict[str, object]:
    """Mean time for a separation delta*size to reach reach*size from points on the attractor, and the slope of that
    time against ln(1/delta); the law says slope * lambda_1 = 1."""
    F = _F(real, over)
    rng = np.random.default_rng(seed)
    times = []
    t_max = 60.0 / max(lam, 1e-6)
    for de in deltas:
        x = np.array(points, float)
        u = rng.standard_normal(x.shape)
        y = x + u * (de * size / np.linalg.norm(u, axis=1, keepdims=True))
        hit = np.full(len(x), np.nan)
        t = 0.0
        while t < t_max and np.isnan(hit).any():
            x, y, t = _rk4(F, x, dt), _rk4(F, y, dt), t + dt
            sep = np.linalg.norm(y - x, axis=1)
            hit[np.isnan(hit) & (sep > reach * size)] = t
        times.append(float(np.nanmean(hit)) if np.isfinite(hit).any() else None)
    ok = [(np.log(1.0 / d), t) for d, t in zip(deltas, times) if t is not None]
    if len(ok) < 3:
        return {"deltas": list(deltas), "times": times, "slope": None, "law": None}
    slope = float(np.polyfit([a for a, _ in ok], [b for _, b in ok], 1)[0])
    return {"deltas": list(deltas), "times": times, "slope": slope, "law": slope * lam}


def chaos(real: Realization, rng, starts: int = 6, transient: float = 300.0, T: float = 1500.0, span: float = 640.0,
          dt: Optional[float] = None, over=None, tol: float = 0.2, folds: float = 10.0) -> Dict[str, object]:
    """Whether the realization is chaotic, its lambda_1, and the predictability law. Chaotic: lambda_1 > 0 at every
    start, with a spread below half its value, the same value over the first half of the run to within a quarter, and
    at least `folds` e-foldings over the run. The last two exclude the separation of neighbouring neutral cycles, which
    grows linearly, so that its finite-time exponent ln(t)/t is positive and halves when the run doubles. The law holds
    if slope * lambda_1 is within tol of 1."""
    dt = dt or min(0.01, real.dt)
    x0 = real.carrier.sample(rng, starts)
    pts = settle(real, x0, transient, dt, over)
    out: Dict[str, object] = {"starts": starts, "escaped": int(starts - len(pts))}
    if not len(pts):
        out.update(chaotic=False, verdict="every start escapes"); return out
    ly = lyapunov(real, pts, T, dt, over=over)
    out["lyapunov"] = ly
    lam = ly.get("median")
    half = None if ly.get("first_half") is None else float(np.median(ly["first_half"]))
    chaotic = bool(lam is not None and lam > 0 and min(ly["lambda1"]) > 0 and ly["spread"] < 0.5 * lam
                   and half is not None and abs(half - lam) < 0.25 * lam and lam * T >= folds)
    if not chaotic:
        out["chaotic"] = False
        out["verdict"] = "not chaotic" if lam is not None else "the exponent could not be measured"
        return out
    sample = [pts[0]]
    F = _F(real, over)
    x = pts[:1].copy()
    every = max(1, int(span / dt) // 64)
    for k in range(int(span / dt)):
        x = _rk4(F, x, dt)
        if k % every == 0:
            sample.append(x[0].copy())
    sample = np.array(sample)
    size = float(np.linalg.norm(sample.max(0) - sample.min(0)))
    out.update(attractor_size=size, chaotic=bool(size > 0))
    if size <= 0:
        out["verdict"] = "not chaotic"
        return out
    # the separation time fluctuates from point to point on the attractor: average over many points along it
    along = sample[:: max(1, len(sample) // 64)][:64]
    pr = predictability(real, along, lam, size, dt, over=over)
    out["predictability"] = pr
    holds = pr["law"] is not None and abs(pr["law"] - 1.0) <= tol
    out["law_holds"] = bool(holds)
    out["verdict"] = (f"chaotic: the initial state is lost at lambda_1 = {lam:.4g} "
                      f"({lam / np.log(2):.3g} bits per unit time); predictability slope x lambda_1 = "
                      f"{pr['law']:.3f}" if pr["law"] is not None else "chaotic; predictability not measured")
    return out


# ------------------------------------------------------------------------------------------------ blow-up
def monomials(real: Realization):
    """Exponent vectors and coefficients of a polynomial drift, per equation; None if the drift is not polynomial."""
    import sympy as sp
    spec = getattr(real, "spec", None)
    if not spec or spec.get("kind") != "equations":
        return None
    from .spec import parse_expression
    variables = list(spec["carrier"]["variables"])
    vs = [sp.Symbol(v, real=True) for v in variables]
    names = {**dict(zip(variables, vs)), **{p: sp.Float(real.params[p]) for p in real.params}}
    eqs = []
    for v in variables:
        e = sp.expand(parse_expression(spec["drift"][v], names))
        try:
            P = sp.Poly(e, *vs)
        except sp.PolynomialError:
            return None
        if any(c.free_symbols for c in P.coeffs()):
            return None
        eqs.append([(np.array(m, float), float(c)) for m, c in P.terms() if abs(float(c)) > 0])
    return {"variables": variables, "equations": eqs}


def balances(mono, max_combinations: int = 200000) -> List[Dict[str, object]]:
    """Dominant balances x_j ~ c_j tau^(-p_j) of a polynomial drift with a growing state. A subset S of the variables
    grows (p_j > 0); in each of their equations one monomial dominates, max_k m_k.p = p_i + 1, and
    p_i c_i = sum of the dominant terms has a real solution with c_i != 0. The other variables have no power growth
    (p_j = 0): their right-hand sides grow at most as tau^(-1), max_k m_k.p <= 1. Below 1 they stay bounded; at 1
    they grow as ln(1/tau), which changes the growing variables only by logarithmic factors (x' = x^3, y' = x^2 - y:
    y ~ (1/2) ln(1/tau)); these are listed as "logarithmic". They enter the growing equations as constants (taken as
    1 in the coefficient equations, which decide only whether a real solution exists)."""
    from scipy.optimize import least_squares
    eqs = mono["equations"]
    d = len(eqs)
    if any(len(e) == 0 for e in eqs):
        return []
    found, seen = [], set()
    for size in range(1, d + 1):
        for S in itertools.combinations(range(d), size):
            choices = [range(len(eqs[i])) for i in S]
            if int(np.prod([len(eqs[i]) for i in S])) > max_combinations:
                continue
            for pick in itertools.product(*choices):
                # unknowns: p_i for i in S; p_j = 0 outside S
                A = np.zeros((size, size))
                for r, (i, k) in enumerate(zip(S, pick)):
                    m = eqs[i][k][0]
                    for c_, j in enumerate(S):
                        A[r, c_] = m[j] - (1.0 if j == i else 0.0)
                try:
                    q = np.linalg.solve(A, np.ones(size))
                except np.linalg.LinAlgError:
                    continue
                if not np.all(np.isfinite(q)) or np.any(q <= 1e-9):
                    continue
                p = np.zeros(d)
                p[list(S)] = q
                if not all(max(m @ p for m, _ in eqs[i]) <= p[i] + 1 + 1e-9 for i in S):
                    continue
                outside = {j: max(m @ p for m, _ in eqs[j]) for j in range(d) if j not in S}
                if any(v > 1 + 1e-9 for v in outside.values()):
                    continue
                log_growing = [int(j) for j, v in outside.items() if v > 1 - 1e-9]
                key = tuple(np.round(p, 6))
                if key in seen:
                    continue
                dom = {i: [(m, a) for m, a in eqs[i] if abs(m @ p - (p[i] + 1)) < 1e-9] for i in S}

                def resid(cS):
                    c = np.ones(d)
                    c[list(S)] = cS
                    return [p[i] * c[i] - sum(a * np.prod(np.abs(c) ** m * np.sign(c) ** m) for m, a in dom[i]) for i in S]

                rng = np.random.default_rng(0)
                best = None
                for _ in range(48):
                    # starts spread over three decades of magnitude: the trivial solution c = 0 attracts small starts
                    c0 = rng.choice([-1.0, 1.0], size) * 10.0 ** rng.uniform(-1, 2, size)
                    sol = least_squares(resid, c0, xtol=1e-12, ftol=1e-12)
                    if np.max(np.abs(sol.fun)) < 1e-8 and np.all(np.abs(sol.x) > 1e-6):
                        best = sol.x
                        break
                if best is None:
                    continue
                seen.add(key)
                c = np.zeros(d)
                c[list(S)] = best
                found.append({"p": p.tolist(), "growing": [int(i) for i in S], "logarithmic": log_growing,
                              "growth": float(p.max()), "n_eff": float(1.0 + 1.0 / p.max()), "c": c.tolist()})
    found.sort(key=lambda b: -b["growth"])
    return found


def escape(real: Realization, x0: np.ndarray, t_max: float = 2000.0, dt: float = 0.01, top: float = 1e8, over=None):
    """Integrate one start with a step limited by the relative change of the state; the escape time (or None) and the
    measured exponent d ln|F| / d ln|x| over amplitudes 1e3 to top."""
    F = _F(real, over)
    x, t, hist = np.array(x0, float)[None], 0.0, []
    while t < t_max:
        f = F(x)
        nx = float(np.linalg.norm(x))
        if float(np.linalg.norm(f)) < 1e-9 * max(1.0, nx):      # settled at an equilibrium: it will not escape
            break
        h = dt if nx < 1 else min(dt, 0.01 * nx / max(float(np.linalg.norm(f)), 1e-300))
        x = _rk4(F, x, h); t += h
        if not np.all(np.isfinite(x)):
            break
        nx = float(np.linalg.norm(x))
        if nx > 1e3:
            hist.append((np.log(nx), np.log(max(float(np.linalg.norm(F(x))), 1e-300))))
        if nx > top:
            break
    escaped = bool(not np.all(np.isfinite(x)) or float(np.linalg.norm(x)) > top)
    n_eff = None
    if len(hist) > 10:
        a = np.array(hist)
        n_eff = float(np.polyfit(a[:, 0], a[:, 1], 1)[0])
    return {"escaped": escaped, "time": t if escaped else None, "n_eff": n_eff}


def blowup(real: Realization, rng, starts: int = 8, over=None, dt: Optional[float] = None,
           tol: float = 0.1) -> Dict[str, object]:
    """Escape of random starts, the measured growth exponent, the escape-time scaling of scaled starts, and the dominant
    balances of the drift that predict them."""
    dt = dt or min(0.01, real.dt)
    x0 = real.carrier.sample(rng, starts)
    runs = [escape(real, x, dt=dt, over=over) for x in x0]
    out: Dict[str, object] = {"escaped": sum(r["escaped"] for r in runs), "starts": starts}
    n_meas = [r["n_eff"] for r in runs if r["escaped"] and r["n_eff"] is not None]
    mono = monomials(real)
    out["balances"] = balances(mono) if mono else None
    if not n_meas:
        out["verdict"] = "no start escapes" if not out["escaped"] else "escapes; growth not measured"
        return out
    n = float(np.median(n_meas))
    out["n_eff"] = n
    if n <= 1.0 + tol:
        out["verdict"] = f"growth without blow-up in finite time (exponent {n:.3g})"
        return out
    p_meas = 1.0 / (n - 1.0)
    out["growth_exponent"] = p_meas
    match = None
    for b in out["balances"] or []:
        if abs(b["n_eff"] - n) <= tol:
            match = b
            break
    out["predicted_by_balance"] = match
    deg = None
    if mono:
        deg = int(max(m.sum() for e in mono["equations"] for m, _ in e))
        out["degree"] = deg
    # escape time of scaled starts. The dominant part of the drift is unchanged by x_j -> s^(p_j) x_j, t -> t / s, so
    # starts scaled with the weights p of the balance escape after a time proportional to 1/s; without a balance, a
    # uniform scaling of the start and the measured exponent give the slope 1 - n_eff instead.
    k = int(np.argmin([r["time"] if r["escaped"] else np.inf for r in runs]))
    base = x0[k] / np.linalg.norm(x0[k])
    if match is not None:
        w = np.clip(np.array(match["p"]), 0.0, None)
        starts_s = [np.sign(base) * np.maximum(np.abs(base), 0.1) * (10.0 * s) ** w for s in (1.0, 2.0, 4.0, 8.0)]
        expected = -1.0
    else:
        starts_s = [1e2 * s * base for s in (1.0, 2.0, 4.0, 8.0)]
        expected = 1.0 - n
    times = [escape(real, x, dt=dt, over=over)["time"] for x in starts_s]
    if all(t is not None for t in times):
        slope = float(np.polyfit(np.log([1, 2, 4, 8]), np.log(times), 1)[0])
        out["escape_time_scaling"] = {"weighted": match is not None, "times": times, "slope": slope, "expected": expected}
    out["verdict"] = (f"blow-up in finite time: |x| ~ (t* - t)^(-{p_meas:.3g})"
                      + (f", predicted by the balance p = {match['growth']:.3g}" if match else ", no balance predicts it"))
    return out


# ------------------------------------------------------------------------------------------------ both
def characterize(real: Realization, rng, starts: int = 4, horizon: float = 100.0, over=None) -> Dict[str, object]:
    """Which of the two, if either, describes a material without a stable state. A first pass integrates a few starts
    together for a fixed horizon: if any escapes, the blow-up test runs; otherwise a short Lyapunov measurement decides
    whether the motion is chaotic, and only then the predictability law is measured. Returns the kind
    ("blow-up", "growth", "chaotic" or None for neither) with the result of the test that decided it."""
    dt = min(0.01, real.dt)
    x0 = real.carrier.sample(rng, starts)
    kept = settle(real, x0, horizon, dt, over)
    if len(kept) < starts:
        b = blowup(real, rng, starts=starts, over=over, dt=dt)
        v = str(b.get("verdict", ""))
        kind = "blow-up" if v.startswith("blow-up") else ("growth" if v.startswith("growth") else None)
        return {"kind": kind, "blowup": b}
    c = chaos(real, rng, starts=3, transient=horizon, T=3.0 * horizon, span=2.0 * horizon, dt=dt, over=over)
    return {"kind": "chaotic" if c.get("chaotic") else None, "chaos": c}
