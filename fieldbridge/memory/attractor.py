"""What a realization without a stable state settles on: chaos, a cycle, a quasi-periodic motion, a drift, convergence
or escape, decided by its two largest Lyapunov exponents and by the recurrence, speed and spectrum of its motion.

The classification is the one the automatic-discovery engine arrived at through its held-out sets (KnowledgeParser,
modules/discovery, corrections 6 to 10 of PROTOCOL.md), ported here so that the record of a realization can carry it.

Exponents. The two largest Lyapunov exponents are measured by the method of Benettin and colleagues: separations
started along random orthonormal directions (along the axes they stay inside one of two uncoupled subsystems),
renormalized by Gram-Schmidt at every step of a fourth-order Runge-Kutta integration at the declared time step, from
STARTS starts that survive a transient; starts that leave every bound are dropped, and the motion is unbounded only
when no start of CANDIDATES stays bounded. The tolerance on an exponent is TOL of the mean contraction rate per
direction, |<tr J>| / n, on the attractor, and at least the error of a finite-time exponent over the run, 25 / T.
Chaotic: lambda_1 above the tolerance over the run and over each half. Converges: lambda_1 below minus the tolerance.
Otherwise the exponents say cycle (lambda_2 below minus the tolerance) or quasi-periodic (both within it), and the
motion of the first start decides.

Motion. A drift along a contracting valley has the exponents of a cycle, and an algebraic approach to rest has
exponents within the tolerance, so the exponents alone are not enough. Converges: the mean speed over the last tenth
of the run is below CONV_FLOOR of the first tenth, or decreases through the tenths of the second half to below
CONV_DROP of its value at the middle, and the second half stays within the range of the first half in every
coordinate, widened by CONV_MARGIN of that range (a motion that slows while it travels, as x ~ t^(1/2), is a drift).
Drift: fewer than REC_MIN upward crossings of the plane through the mean of the most varying coordinate. Cycle or
quasi-periodic: the record is extended to SPEC_PERIODS returns and the sharp lines of the spectrum of one generic
coordinate are compared; two lines are independent when their ratio lies farther than the uncertainty of SPEC_BINS
bins from every rational with denominator up to SPEC_DENOM. Two independent frequencies with lambda_1 within the
tolerance over each half make a quasi-periodic motion; one frequency a cycle; no sharp line is undetermined.

Numerical statements about the supplied equations. A chaotic verdict needs lambda_1 above its tolerance on every start;
the laws of chaotic loss and of blow-up are in chaos_blowup.
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np

from . import analysis as an
from .spec import Realization

STARTS, TRANSIENT, TOL, STEPS = 4, 5000, 0.01, 400_000
CANDIDATES = 16
REC_MIN, REC_PERIOD, REC_TOL = 20, 8, 1e-3
CONV_FLOOR, CONV_DROP, CONV_MARGIN = 1e-3, 0.9, 0.01
SPEC_PERIODS, SPEC_STEPS, SPEC_FLOOR, SPEC_PEAKS, SPEC_BINS, SPEC_DENOM = 1000, 1_000_000, 1e-3, 8, 0.25, 12


def _rk4(real: Realization, X: np.ndarray, dt: float) -> np.ndarray:
    k1 = real.F(X)
    k2 = real.F(real.carrier.wrap(X + 0.5 * dt * k1))
    k3 = real.F(real.carrier.wrap(X + 0.5 * dt * k2))
    k4 = real.F(real.carrier.wrap(X + dt * k3))
    return (k1 + 2 * k2 + 2 * k3 + k4) / 6.0


def attractor(real: Realization, rng, steps_cap: int = STEPS, spec_steps: int = SPEC_STEPS) -> Dict[str, object]:
    """The kind of motion of a realization without a stable state, with the exponents and the tests that decided it.
    steps_cap bounds the exponent run and spec_steps the record of the spectrum (the defaults are those of the
    engine; a slow escape that outlasts the run reads as a drift, which the blow-up law of chaos_blowup resolves)."""
    n = real.carrier.dim
    dt = an.safe_dt(real)
    step = lambda X: real.carrier.wrap(X + dt * _rk4(real, X, dt))  # noqa: E731
    period = getattr(real.carrier, "period", None) if real.carrier.kind == "torus" else None
    bound = 1e6 * max(1.0, an.base_length(real))

    def bounded(X):
        with np.errstate(all="ignore"):
            return np.all(np.isfinite(X), axis=-1) & (np.max(np.abs(np.nan_to_num(X, nan=np.inf)), axis=-1) <= bound)

    def advance(X, k):
        for _ in range(k):
            if X.shape[0] == 0:
                break
            X = step(X)
            X = X[bounded(X)]
        return X

    X = advance(real.carrier.sample(rng, CANDIDATES), TRANSIENT)
    escaped = CANDIDATES - X.shape[0]
    if X.shape[0] == 0:
        return {"kind": "unbounded", "escaped": escaped}
    X = X[:STARTS]
    traces = []
    probe = X
    for _ in range(20):
        probe = advance(probe, 50)
        if probe.shape[0] == 0:
            return {"kind": "unbounded", "escaped": CANDIDATES}
        traces += [np.trace(an.jacobian(real, x)) for x in probe]
    rate = abs(float(np.mean(traces))) / n
    # a second transient of twenty contraction times, so that slow approaches end before the exponents are measured
    X = advance(probe, int(min(steps_cap // 4, 20.0 / max(rate, 1e-12) / dt)))
    if X.shape[0] == 0:
        return {"kind": "unbounded", "escaped": CANDIDATES}
    tol = TOL * max(rate, 1e-12)
    steps = int(min(steps_cap, 50.0 / tol / dt))
    capped = 50.0 / tol / dt > steps_cap
    tol = max(tol, 25.0 / (steps * dt))
    delta = 1e-7 * max(1.0, an.base_length(real))
    m = min(2, n)
    S = X.shape[0]
    E = np.stack([np.linalg.qr(rng.standard_normal((n, m)))[0].T for _ in range(S)])
    logs = np.zeros((2, S, m))
    path = []  # the first measured start, for the tests of the motion
    for k in range(steps):
        B = np.concatenate([X[:, None, :], X[:, None, :] + delta * E], axis=1).reshape(-1, n)
        B = step(B).reshape(X.shape[0], m + 1, n)
        alive = np.all(bounded(B), axis=1)
        if not alive.all():  # a start that leaves every bound is dropped with its exponents
            B, E, logs = B[alive], E[alive], logs[:, alive]
            path = []
            if B.shape[0] == 0:
                return {"kind": "unbounded", "escaped": CANDIDATES}
        X, D = B[:, 0], B[:, 1:] - B[:, :1]
        path.append(X[0].copy())
        if period is not None:
            D = (D + period / 2) % period - period / 2
        half = int(2 * k >= steps)
        for i in range(m):  # Gram-Schmidt on the separations
            for j in range(i):
                D[:, i] -= np.sum(D[:, i] * E[:, j], axis=1, keepdims=True) * E[:, j]
            norm = np.linalg.norm(D[:, i], axis=1)
            logs[half, :, i] += np.log(norm / delta)
            E[:, i] = D[:, i] / norm[:, None]
    T = steps * dt
    full = np.median(logs.sum(axis=0) / T, axis=0)
    halves = np.median(logs / (T / 2), axis=1)
    l1, l2 = float(full[0]), (float(full[1]) if m > 1 else None)
    if l1 > tol and np.all(halves[:, 0] > tol):
        kind = "chaotic"
    elif l1 < -tol:
        kind = "converges"
    elif abs(l1) <= tol and (l2 is None or l2 < -tol):
        kind = "cycle"
    elif abs(l1) <= tol and l2 is not None and abs(l2) <= tol:
        kind = "quasi-periodic"
    else:
        kind = "undetermined"
    out: Dict[str, object] = {"kind": kind, "lyapunov": [l1, l2], "halves": halves[:, 0].tolist(), "tolerance": tol,
                              "contraction_rate": rate, "time": T, "capped": capped, "starts": int(logs.shape[1]),
                              "escaped": escaped}
    if kind in ("cycle", "quasi-periodic") and len(path) > 2:
        out["exponents_kind"] = kind
        out.update(_motion(real, np.asarray(path), step, dt, halves[:, 0], tol, rng, spec_steps))
    return out


def _motion(real: Realization, P: np.ndarray, step, dt: float, halves, tol: float, rng,
            spec_steps: int = SPEC_STEPS) -> Dict[str, object]:
    """Converges, drift, cycle, quasi-periodic or undetermined, for a motion whose largest exponent lies within the
    tolerance (see the module docstring)."""
    if len(P) < 10 * REC_MIN:
        return {"kind": "undetermined", "record": "too short"}
    with np.errstate(all="ignore"):
        speed = np.linalg.norm(real.F(P), axis=1)
    tenths = np.array([float(np.mean(s)) for s in np.array_split(speed, 10)])
    out: Dict[str, object] = {"speed_tenths": tenths.tolist()}
    slows = bool(tenths[-1] < CONV_FLOOR * tenths[0]
                 or (np.all(np.diff(tenths[5:]) < 0) and tenths[-1] < CONV_DROP * tenths[5]))
    first, second = P[: len(P) // 2], P[len(P) // 2:]
    lo, hi = first.min(axis=0), first.max(axis=0)
    margin = CONV_MARGIN * (hi - lo)
    out["in_place"] = bool(np.all(second.min(axis=0) >= lo - margin) and np.all(second.max(axis=0) <= hi + margin))
    if slows and out["in_place"]:
        out["kind"] = "converges"
        return out
    out["recurrence"] = recurrence(P)
    returns = out["recurrence"]["crossings"]
    if returns < REC_MIN:
        out["kind"] = "drift"
        return out
    extra = int(min(spec_steps - len(P), len(P) * (SPEC_PERIODS / returns - 1)))
    if extra > 0:  # the record of the first start, continued until it holds SPEC_PERIODS returns
        X, more, k = P[-1:].copy(), np.empty((extra, P.shape[1])), 0
        while k < extra:
            X = step(X)
            if not np.all(np.isfinite(X)):
                break
            more[k], k = X[0], k + 1
        P = np.concatenate([P, more[:k]])
    out["spectrum"] = frequencies(P, dt, rng)
    independent = out["spectrum"]["independent"]
    within = bool(np.all(np.abs(halves) <= tol))
    out["kind"] = ("quasi-periodic" if independent == 2 and within else "cycle" if independent == 1 else "undetermined")
    return out


def frequencies(P: np.ndarray, dt: float, rng) -> Dict[str, object]:
    """Sharp lines in the spectrum of one generic coordinate (the standardized coordinates projected on a random
    direction, since no single coordinate of two uncoupled oscillators carries both frequencies): a line is a local
    maximum of the Hann-windowed periodogram above SPEC_FLOOR of the strongest, with at least REC_MIN periods in the
    record and the power two bins away below a tenth of its own, its frequency interpolated from the logarithm of the
    power. The strongest line b and each of the next lines g: g is commensurate with b when |g/b - p/q| lies within
    SPEC_BINS (1 + g/b)/(b T) for some q up to SPEC_DENOM; the first g that is not makes two independent frequencies."""
    sd = P.std(axis=0)
    live = sd > 1e-12 * max(float(sd.max()), 1e-300)
    w = rng.standard_normal(P.shape[1]) * live
    s = ((P - P.mean(axis=0)) / np.where(live, sd, 1.0)) @ (w / max(float(np.linalg.norm(w)), 1e-300))
    T = s.size * dt
    W = np.abs(np.fft.rfft((s - s.mean()) * np.hanning(s.size))) ** 2
    if W.size < REC_MIN + 3 or not np.any(W[REC_MIN:] > 0):
        return {"independent": 0, "lines": [], "time": T}
    i = np.arange(REC_MIN, W.size - 2)
    peak = (W[i] > W[i - 1]) & (W[i] >= W[i + 1]) & (W[i] > SPEC_FLOOR * W[REC_MIN:].max())
    sharp = (W[i - 2] < 0.1 * W[i]) & (W[i + 2] < 0.1 * W[i])
    i = i[peak & sharp]
    i = i[np.argsort(-W[i])][:SPEC_PEAKS]
    if i.size == 0:
        return {"independent": 0, "lines": [], "time": T}
    a, b0, c = (np.log(np.maximum(W[i + k], 1e-300)) for k in (-1, 0, 1))
    lines = (i + 0.5 * (a - c) / (a - 2 * b0 + c)) / T
    base, out = float(lines[0]), {"independent": 1, "lines": lines.tolist(), "time": T}
    for g in lines[1:]:
        r = float(g) / base
        tol_r = SPEC_BINS * (1.0 + r) / (base * T)
        if all(abs(q * r - round(q * r)) > q * tol_r for q in range(1, SPEC_DENOM + 1)):
            out.update(independent=2, ratio=r, ratio_tolerance=tol_r)
            break
    return out


def recurrence(P: np.ndarray) -> Dict[str, object]:
    """Upward crossings of the plane through the mean of the most varying coordinate. Fewer than REC_MIN crossings:
    the motion does not return (drift). Otherwise the returns of the second half are split in two blocks; the period
    n up to REC_PERIOD that brings a return closest to the n-th next one gives the distance d relative to the size of
    the motion, and a cycle has d within REC_TOL in the last block or halved from the first block to the last."""
    k = int(np.argmax(P.std(axis=0)))
    level = float(P[:, k].mean())
    up = np.flatnonzero((P[:-1, k] < level) & (P[1:, k] >= level))
    if up.size < REC_MIN:
        return {"kind": "drift", "crossings": int(up.size)}
    size = max(float(np.max(P.max(axis=0) - P.min(axis=0))), 1e-300)
    R = P[up + 1][up.size // 2:]
    q = len(R) // 2

    def closest(block):
        best, per = np.inf, None
        for p in range(1, min(REC_PERIOD, len(block) - 2) + 1):
            d = float(np.mean(np.linalg.norm(block[p:] - block[:-p], axis=1))) / size
            if d < best:
                best, per = d, p
        return best, per

    d_first, _ = closest(R[:q])
    d_last, per = closest(R[q:])
    cycle = d_last <= REC_TOL or d_last < 0.5 * d_first
    return {"kind": "cycle" if cycle else "quasi-periodic", "crossings": int(up.size), "distance": [d_first, d_last],
            "period": per}
