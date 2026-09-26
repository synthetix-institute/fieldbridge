"""Memory in the phase of an oscillator.

A realization without a stable state can still remember. On a limit cycle, time-translation symmetry makes the
phase a flat direction: a pulse shifts the phase, the shift persists because nothing restores it, and noise loses
it by phase diffusion, Law 2. A bounded write that shifts the phase at a rate proportional to its strength then
obeys the flat-direction form of the relation between retention and writing times: their ratio grows linearly with
the write strength.

Phases are read from the times at which the first coordinate crosses a section upward (its mean over the cycle).
The k-th crossing of a trajectory, compared with k periods of the unperturbed cycle, gives its phase lag.
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence

import numpy as np

from . import analysis as an
from .identity import Realization


def _crossings(real: Realization, q: np.ndarray, T: float, D: float, rng, section: float, dt: float,
               over: Optional[Dict] = None, max_cross: int = 400, forcing=None, band: float = 0.0):
    """Integrate an ensemble and return the upward crossing times of q[:, 0] through `section` (NaN-padded).
    A crossing counts only after the coordinate has been below section - band since the previous crossing, so
    that noise near the section does not add crossings."""
    over = dict(over or {})
    n = q.shape[0]
    times = np.full((n, max_cross), np.nan)
    count = np.zeros(n, int)
    armed = q[:, 0] < section - band
    amp = np.sqrt(2.0 * D * dt)
    prev = q[:, 0].copy()
    steps = int(round(T / dt))
    for k in range(steps):
        t = k * dt
        if forcing is not None:
            over.update(forcing(t))
        q = real.carrier.wrap(q + dt * real.F(q, **over) + amp * rng.standard_normal(q.shape))
        cur = q[:, 0]
        armed |= cur < section - band
        up = armed & (prev < section) & (cur >= section) & (count < max_cross)
        if np.any(up):
            idx = np.flatnonzero(up)
            frac = (section - prev[idx]) / np.maximum(cur[idx] - prev[idx], 1e-300)
            times[idx, count[idx]] = t + frac * dt
            count[idx] += 1
            armed[idx] = False
        prev = cur.copy()
    return q, times


def limit_cycle(real: Realization, rng, t_transient: float = 300.0, t_record: float = 200.0) -> Dict[str, object]:
    """Period, section and a point on the cycle, from one noise-free trajectory."""
    dt = an.safe_dt(real)
    q0 = real.carrier.sample(rng, 1)
    q, _, _ = an.integrate(real, q0, t_transient, 0.0, rng)
    _, t, rec = an.integrate(real, q, t_record, 0.0, rng, record=lambda x: x[0].copy(), n_record=int(t_record / dt))
    rec = np.asarray(rec)
    section = float(rec[:, 0].mean())
    _, times = _crossings(real, q, t_record, 0.0, rng, section, dt)
    ts = times[0][np.isfinite(times[0])]
    if ts.size < 3:
        raise ValueError("no sustained oscillation: the realization settles, and its memory is not a phase")
    period = float(np.mean(np.diff(ts)))
    # a point on the cycle exactly at the section crossing
    q_at, _ = _advance_to_crossing(real, q, section, dt)
    return {"period": period, "section": section, "state": q_at[0].tolist(), "amplitude": float(np.ptp(rec[:, 0]))}


def _advance_to_crossing(real, q, section, dt):
    prev = q[:, 0].copy()
    for k in range(10 ** 6):
        q = real.carrier.wrap(q + dt * real.F(q))
        if prev[0] < section <= q[0, 0]:
            return q, k * dt
        prev = q[:, 0].copy()
    raise RuntimeError("no crossing")


def phase_lag(real: Realization, q: np.ndarray, cycle: Dict, rng, n_periods: int = 6) -> np.ndarray:
    """Asymptotic phase lag of states q relative to the reference crossing, in radians (noise-free)."""
    dt = an.safe_dt(real)
    T = cycle["period"]
    _, times = _crossings(real, q, n_periods * T, 0.0, rng, cycle["section"], dt)
    last = np.array([row[np.isfinite(row)][-1] for row in times])
    k = np.array([np.isfinite(row).sum() for row in times])
    return 2 * np.pi * ((last - k * T) / T)


def response_curve(real: Realization, cycle: Dict, rng, kick: Sequence[float], n: int = 24) -> Dict[str, object]:
    """Phase shift left by an instantaneous kick applied at n phases of the cycle (the phase response curve)."""
    dt = an.safe_dt(real)
    T = cycle["period"]
    q = np.asarray(cycle["state"], float)[None]
    states, phases = [], []
    for j in range(n):
        states.append(q[0].copy())
        phases.append(2 * np.pi * j / n)
        q, _, _ = an.integrate(real, q, T / n, 0.0, rng, dt=dt)
    S = np.array(states)
    base = phase_lag(real, S, cycle, rng)
    kicked = phase_lag(real, real.carrier.wrap(S + np.asarray(kick, float)[None]), cycle, rng)
    shift = np.angle(np.exp(1j * (kicked - base)))
    return {"phase": phases, "shift": shift.tolist()}


def hold(real: Realization, cycle: Dict, rng, shift: float, D: float, n_traj: int = 300, n_periods: int = 60,
         phase_of_write: float = 0.0) -> Dict[str, object]:
    """Two ensembles, one written (phase advanced by `shift`) and one not, under noise D: the mean lag stays,
    its variance grows linearly (Law 2), and the separation d' = shift / std decays as t^{-1/2}."""
    dt = an.safe_dt(real)
    T = cycle["period"]
    q = np.repeat(np.asarray(cycle["state"], float)[None], n_traj, axis=0)
    q_written, _, _ = an.integrate(real, q[:1], (shift / (2 * np.pi)) * T, 0.0, rng, dt=dt)
    qa = q.copy()
    qb = np.repeat(q_written, n_traj, axis=0)
    band = 0.2 * cycle["amplitude"]
    _, ta = _crossings(real, qa, n_periods * T, D, rng, cycle["section"], dt, band=band)
    _, tb = _crossings(real, qb, n_periods * T, D, rng, cycle["section"], dt, band=band)
    K = min(int(np.min(np.isfinite(ta).sum(1))), int(np.min(np.isfinite(tb).sum(1))))
    ka = np.arange(1, K + 1)
    lag_a = 2 * np.pi * (ta[:, :K] - ka[None] * T) / T
    lag_b = 2 * np.pi * (tb[:, :K] - ka[None] * T) / T
    # the written ensemble crosses earlier by `shift`
    sep = np.mean(lag_a, 0) - np.mean(lag_b, 0)
    var = 0.5 * (np.var(lag_a, 0) + np.var(lag_b, 0))
    t = ka * T
    slope = float(np.polyfit(t, var, 1)[0])
    return {"time": t.tolist(), "separation": sep.tolist(), "variance": var.tolist(),
            "phase_diffusion": slope / 2.0, "d_prime": (sep / np.sqrt(np.maximum(var, 1e-300))).tolist()}


def lock(real: Realization, cycle: Dict, rng, strengths: Sequence[float], target: float, D: float,
         n_traj: int = 200) -> Dict[str, object]:
    """Rewrite the phase by a weak sustained change of the control (a frequency shift) of relative strength h, and
    compare with the retention time of a shift of the same size. Along a flat direction T_hold / T_write grows linearly
    in h."""
    if not real.control:
        raise ValueError("a control parameter is needed to shift the frequency")
    dt = an.safe_dt(real)
    T = cycle["period"]
    v0 = float(real.params[real.control])
    q = np.repeat(np.asarray(cycle["state"], float)[None], n_traj, axis=0)
    out = []
    hold_run = hold(real, cycle, rng, target, D, n_traj=n_traj, n_periods=40)
    Dphi = hold_run["phase_diffusion"]
    T_hold = target ** 2 / (2.0 * Dphi)
    for h in strengths:
        over = {real.control: v0 * (1.0 + h)}
        # period of the shifted cycle, noise-free
        _, times = _crossings(real, q[:1], 12 * T, 0.0, rng, cycle["section"], dt, over=over)
        ts = times[0][np.isfinite(times[0])]
        T_h = float(np.mean(np.diff(ts[2:])))
        rate = abs(2 * np.pi / T_h - 2 * np.pi / T)
        T_write = target / rate
        out.append({"h": float(h), "period": T_h, "phase_rate": rate, "T_write": T_write,
                    "ratio": T_hold / T_write})
    hs = np.array([o["h"] for o in out])
    rs = np.array([o["ratio"] for o in out])
    slope = float(np.polyfit(np.log(hs), np.log(rs), 1)[0])
    return {"target_shift": target, "phase_diffusion": Dphi, "T_hold": T_hold, "writes": out,
            "log_log_slope": slope}
