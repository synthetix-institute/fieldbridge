#!/usr/bin/env python3
"""Finite populations of agents against the mean-field equations of examples/memory/populations.

Each model is simulated as individual events in a population of size N, with no reference to its mean-field equation:

  kuramoto   N phase oscillators with Lorentzian natural frequencies, coupled through their mean field
  cycle      N agents playing three strategies in cyclic dominance: pairwise comparison and mutation
  error      N replicating sequences: a parent chosen by fitness, copied with fidelity q, replacing a random agent
  budworm    individuals with birth, crowding and predation, in a habitat of size N

and the results are compared with what the equations say at the same parameters: the place of a transition, the state
beyond it, the relaxation rate near it, and how the time to leave a state grows with N. For the two birth-death models
the mean time to leave a state is also known exactly from the rates, which tests the simulation itself. The part
"closed" gathers what follows from the equations and the exact chain alone, with repeated sweeps at small N to test it.

Usage: python3 -B compare.py [--parts kuramoto cycle error budworm closed fold_passage actions] [--quick] [--out results.json]
Needs numpy and scipy; numba makes the event loops about a hundred times faster (pip install numba). On four cores the
full run takes about four hours, nearly all of it the budworm sweeps at N = 100000; --quick leaves those out.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.signal import welch
from scipy.special import logsumexp

try:
    from numba import njit
except ImportError:  # the same code runs without it, more slowly
    def njit(f):
        return f


# ------------------------------------------------------------------------------------------------ Kuramoto
def lorentzian(n, delta, rng):
    """Natural frequencies at the quantiles of a Lorentzian: no sample mean to drift with, and the tails are kept."""
    q = (np.arange(n) + 0.5) / n
    return delta * np.tan(np.pi * (q - 0.5))


def kuramoto_run(n, K, delta, T, dt, rng, theta0=None, noise=0.0, every=10):
    w = lorentzian(n, delta, rng)
    th = rng.uniform(0, 2 * np.pi, n) if theta0 is None else np.full(n, float(theta0))
    steps, out = int(T / dt), []
    amp = np.sqrt(2 * noise * dt)
    for k in range(steps):
        z = np.exp(1j * th).mean()
        if k % every == 0:
            out.append(z)
        # midpoint step
        mid = th + 0.5 * dt * (w + K * np.imag(z * np.exp(-1j * th)))
        zm = np.exp(1j * mid).mean()
        th = th + dt * (w + K * np.imag(zm * np.exp(-1j * mid)))
        if noise:
            th = th + amp * rng.standard_normal(n)
    return np.array(out), every * dt


def part_kuramoto(rng, quick):
    delta = 0.5
    Ns = (100, 1000, 10000)
    Ks = np.array([0.4, 0.7, 0.9, 1.0, 1.1, 1.3, 1.6, 2.0, 2.5, 3.0])
    res = {"delta": delta, "K": Ks.tolist(), "theory_r": [float(np.sqrt(max(0.0, 1 - 2 * delta / K))) for K in Ks],
           "N": list(Ns), "r": {}, "r_sd": {}}
    for n in Ns:
        rs, sds = [], []
        for K in Ks:
            z, _ = kuramoto_run(n, K, delta, 200.0 if n > 5000 else 400.0, 0.02, rng)
            r = np.abs(z[len(z) // 2:])
            rs.append(float(r.mean()))
            sds.append(float(r.std()))
        res["r"][str(n)] = rs
        res["r_sd"][str(n)] = sds
    # below the transition the mean field is a finite-size fluctuation of size N^(-1/2)
    i = 0
    res["fluctuation_exponent"] = float(np.polyfit(np.log(Ns), np.log([res["r"][str(n)][i] for n in Ns]), 1)[0])
    # relaxation of a prepared coherent state: all phases equal lies on the Ott-Antonsen manifold, where r obeys the
    # closed equation exactly; the rate near the transition goes to zero as |K/2 - Delta|
    n = 40000
    rel = []
    for K in (0.2, 0.5, 0.7, 0.8, 0.9):
        z, h = kuramoto_run(n, K, delta, 60.0, 0.02, rng, theta0=0.0)
        r, t = np.abs(z), np.arange(len(z)) * h
        ok = (r < 0.25) & (r > 0.03)
        rate = float(-np.polyfit(t[ok], np.log(r[ok]), 1)[0]) if ok.sum() > 5 else None
        rel.append({"K": K, "rate": rate, "theory": delta - K / 2})
    res["relaxation"] = rel
    # above the transition the collective phase has no restoring force: with noise on each oscillator it diffuses
    # at a rate proportional to 1/N, whereas the modulus stays put
    dif = []
    for n in (100, 300, 1000, 3000):
        slopes = []
        for _ in range(12 if not quick else 6):
            z, h = kuramoto_run(n, 2.0, delta, 300.0, 0.02, rng, noise=0.02)
            z = z[len(z) // 3:]
            psi = np.unwrap(np.angle(z))
            t = np.arange(len(psi)) * h
            psi = psi - np.polyval(np.polyfit(t, psi, 1), t)
            lags = np.arange(1, len(psi) // 4, 25)
            msd = [np.mean((psi[l:] - psi[:-l]) ** 2) for l in lags]
            slopes.append(np.polyfit(lags[:6] * h, msd[:6], 1)[0] / 2)
        dif.append({"N": n, "D": float(np.mean(slopes)), "D_se": float(np.std(slopes) / np.sqrt(len(slopes)))})
    res["phase_diffusion"] = dif
    res["phase_diffusion_exponent"] = float(np.polyfit(np.log([d["N"] for d in dif]), np.log([d["D"] for d in dif]), 1)[0])
    return res


# ------------------------------------------------------------------------------------------------ cyclic dominance
@njit
def cycle_events(n0, N, a, b, mu, T, every, seed):
    """Well-mixed agents, one event at a time. With probability m an agent takes a uniformly drawn strategy; otherwise
    agent j copies agent i with probability (1 + (pi_i - pi_j)/P)/2. In the limit N -> infinity the frequencies obey
    the replicator equation with mutation rate mu in time t = events * (1 - m) / (P N), with P = a + b."""
    np.random.seed(seed)
    P = a + b
    m = mu / (P + mu)
    n = n0.copy()
    dt_event = (1 - m) / (P * N)
    samples = int(T / every)
    out = np.zeros((samples, 3))
    t, next_t, k = 0.0, 0.0, 0
    while k < samples:
        if t >= next_t:
            out[k] = n / N
            k += 1
            next_t += every
        t += dt_event
        u = np.random.random()
        v = int(np.random.random() * N)
        kind_j = 0 if v < n[0] else (1 if v < n[0] + n[1] else 2)
        if u < m:
            new = int(np.random.random() * 3)
            n[kind_j] -= 1
            n[new] += 1
            continue
        v = int(np.random.random() * N)
        kind_i = 0 if v < n[0] else (1 if v < n[0] + n[1] else 2)
        if kind_i == kind_j:
            continue
        x0, x1, x2 = n[0] / N, n[1] / N, n[2] / N
        pay = (a * x1 - b * x2, a * x2 - b * x0, a * x0 - b * x1)
        if np.random.random() < 0.5 * (1 + (pay[kind_i] - pay[kind_j]) / P):
            n[kind_j] -= 1
            n[kind_i] += 1
    return out


def cycle_meanfield(a, b, mu, T=4000.0):
    def f(t, s):
        x, y = s
        z = 1 - x - y
        phi = (a - b) * (x * y + y * z + z * x)
        return [x * (a * y - b * z - phi) + mu * (1 / 3 - x), y * (a * z - b * x - phi) + mu * (1 / 3 - y)]
    sol = solve_ivp(f, (0, T), [0.4, 0.3], method="LSODA", rtol=1e-8, atol=1e-10, t_eval=np.arange(T / 2, T, 0.5))
    return float(np.sqrt(np.mean((sol.y[0] - 1 / 3) ** 2))), spectral_peak(sol.y[0] - 1 / 3, 0.5)


def spectral_peak(d, h):
    """Angular frequency of the largest peak of an averaged spectrum; a single periodogram scatters over a broad peak."""
    w, p = welch(d, 1 / h, nperseg=len(d) // 4, detrend="constant")
    return float(2 * np.pi * w[1:][np.argmax(p[1:])])


def part_cycle(rng, quick):
    a, b = 1.0, 1.5
    mu_h = (b - a) / 6
    omega = np.sqrt(3) * (a + b) / 6
    Ns = (1000, 10000) if quick else (1000, 10000, 100000)
    mus = np.array([0.02, 0.04, 0.06, 0.07, 0.08, 0.0833, 0.09, 0.10, 0.12, 0.16, 0.24])
    T = 300.0 if quick else 1500.0
    res = {"a": a, "b": b, "mu_h": mu_h, "omega": omega, "mu": mus.tolist(), "N": list(Ns), "rms": {}, "peak": {}}
    res["meanfield_rms"], res["meanfield_peak"] = map(list, zip(*(cycle_meanfield(a, b, m) for m in mus)))
    res["meanfield_peak"] = [p if v > 1e-3 else None for v, p in zip(res["meanfield_rms"], res["meanfield_peak"])]
    for N in Ns:
        rms, peak = [], []
        for m in mus:
            n0 = np.array([N // 3, N // 3, N - 2 * (N // 3)], dtype=np.int64)
            tt = T
            x = cycle_events(n0, N, a, b, float(m), tt, 0.5, int(rng.integers(2 ** 31)))[len(np.arange(0, tt, 0.5)) // 5:, 0]
            d = x - 1 / 3
            rms.append(float(np.sqrt(np.mean(d ** 2))))
            peak.append(spectral_peak(d, 0.5))
        res["rms"][str(N)] = rms
        res["peak"][str(N)] = peak
    # above the transition the state at 1/3 is stable: finite N sustains cycles at the frequency of the eigenvalues.
    # Copying and mutation events give the noise intensity B = (2P + 4 mu) / (9 N) on each frequency; it is isotropic in
    # the plane of the simplex, so that it commutes with the rotation and the variance is B / (2 (mu - mu_h))
    P = a + b
    res["linear_noise_variance_times_N"] = [float((P + 2 * m) / (9 * (m - mu_h))) if m > mu_h else None for m in mus]
    i = list(mus).index(0.16)
    res["quasi_cycle_exponent"] = float(np.polyfit(np.log(Ns), np.log([res["rms"][str(N)][i] for N in Ns]), 1)[0])
    res["quasi_cycle_peak_at_mu_0.16"] = {str(N): res["peak"][str(N)][i] for N in Ns}
    return res


# ------------------------------------------------------------------------------------------------ birth-death chains
def exit_time(birth, death, start, stop, N_max=None):
    """Exact mean first-passage time of a birth-death chain from start to stop, with a reflecting state at the far
    side (0 is reflecting when the chain is run upward; rates are per unit time)."""
    if stop > start:   # upward: T(n -> n+1) = (1 + d_n T(n-1 -> n)) / b_n, with reflection at the lowest state
        log_tau, total = None, 0.0
        step = 0.0
        low = 1
        for n in range(low, stop):
            step = (1.0 + death(n) * step) / birth(n)
            if n >= start:
                total += step
        return total
    top = N_max
    step, total = 0.0, 0.0
    for n in range(top, stop, -1):
        step = (1.0 + birth(n) * step) / death(n)
        if n <= start:
            total += step
    return total


def extinction_time(birth, death, N, start):
    """Natural logarithm of the exact mean time to reach 0 from start, for a chain absorbed at 0 and bounded at N. The
    differences D_j = tau_j - tau_(j-1) obey d_j D_j = 1 + b_j D_(j+1), with b_N = 0."""
    D, total = 0.0, 0.0
    for j in range(N, 0, -1):
        D = (1.0 + birth(j) * D) / death(j)
        if j <= start:
            total += D
    return float(np.log(total))


@njit
def error_events(n0, N, f0, q, T, every, seed):
    np.random.seed(seed)
    samples = int(T / every)
    out = np.zeros(samples)
    n = n0
    t, next_t, k = 0.0, 0.0, 0
    while k < samples:
        if t >= next_t:
            out[k] = n / N
            k += 1
            next_t += every
        t += 1.0 / N
        W = f0 * n + (N - n)
        parent_master = np.random.random() < f0 * n / W
        child_master = parent_master and np.random.random() < q
        dead_master = np.random.random() < n / N
        n += int(child_master) - int(dead_master)
    return out


@njit
def error_lifetimes(n0, N, f0, q, reps, tmax, seed):
    np.random.seed(seed)
    out = np.zeros(reps)
    for r in range(reps):
        n, t = n0, 0.0
        while n > 0 and t < tmax:
            t += 1.0 / N
            W = f0 * n + (N - n)
            pm = np.random.random() < f0 * n / W
            cm = pm and np.random.random() < q
            dm = np.random.random() < n / N
            n += int(cm) - int(dm)
        out[r] = t
    return out


def part_error(rng, quick):
    f0 = 2.0
    qc = 1 / f0
    qs = np.array([0.40, 0.44, 0.47, 0.49, 0.50, 0.51, 0.53, 0.56, 0.60, 0.70, 0.80, 0.90])
    Ns = (100, 1000, 10000) if quick else (100, 1000, 10000, 100000)
    res = {"f0": f0, "q_c": qc, "q": qs.tolist(), "theory_x": [max(0.0, (f0 * q - 1) / (f0 - 1)) for q in qs], "N": list(Ns),
           "x": {}}
    for N in Ns:
        xs = []
        for q in qs:
            # the mean over runs started at the master sequence, read after the slowest relaxation near q_c
            T = 400.0 if N <= 10000 else 150.0
            vals = [error_events(N, N, f0, float(q), T, 1.0, int(rng.integers(2 ** 31)))[-int(T / 4):].mean()
                    for _ in range(4)]
            xs.append(float(np.mean(vals)))
        res["x"][str(N)] = xs
    q = 0.8
    birth = lambda N: (lambda n: N * (f0 * n * q / (f0 * n + N - n)) * (1 - n / N))
    death = lambda N: (lambda n: N * (n / N) * ((f0 * n * (1 - q) + (N - n)) / (f0 * n + N - n)))
    # the lifetime of the master sequence: exact for the chain, simulated where it is short enough
    x_star = (f0 * q - 1) / (f0 - 1)
    life = []
    for N in (10, 15, 20, 25, 30, 40, 60, 100, 200):
        lt = extinction_time(birth(N), death(N), N, int(round(x_star * N)))
        row = {"N": N, "ln_tau_exact": lt}
        if lt < np.log(3e4):
            times = error_lifetimes(int(round(x_star * N)), N, f0, q, 300, 1e7, int(rng.integers(2 ** 31)))
            row.update(ln_tau_sim=float(np.log(times.mean())), ln_tau_sim_se=float(times.std() / times.mean() / np.sqrt(len(times))))
        life.append(row)
    res["lifetime"] = life
    # large-deviation action: S = int_0^x* ln(b/d) dx with b, d the rates per agent; ln tau = N S + o(N)
    def ratio(x):
        return (f0 * q * (1 - x)) / (f0 * x * (1 - q) + 1 - x)
    xs = np.linspace(1e-6, x_star, 20001)
    S = float(np.trapezoid(np.log(ratio(xs)), xs) if hasattr(np, "trapezoid") else np.trapz(np.log(ratio(xs)), xs))
    res["action"] = S
    big = [r for r in life if r["N"] >= 60]
    res["action_from_exact"] = float(np.polyfit([r["N"] for r in big], [r["ln_tau_exact"] for r in big], 1)[0])
    # the same model with constant noise, as a specification has it: barrier V / D with D = noise intensity at the state
    drift = lambda x: x * (f0 * q - 1 - (f0 - 1) * x) / (1 + (f0 - 1) * x)
    grid = np.linspace(0, x_star, 20001)
    dV = float(np.trapezoid(drift(grid), grid) if hasattr(np, "trapezoid") else np.trapz(drift(grid), grid))
    b_star = (f0 * q * (1 - x_star)) * x_star / (1 + (f0 - 1) * x_star)
    res["constant_noise_action"] = dV / b_star     # D = (b + d) / (2N) with b = d at the state
    return res


# ------------------------------------------------------------------------------------------------ budworm
@njit
def budworm_events(n0, N, r, k, T, every, seed, k_end):
    """Gillespie events: birth r u, crowding r u^2/k and predation u^2/(1+u^2), per unit habitat, with u = n/N. The
    capacity moves linearly from k to k_end over the run."""
    np.random.seed(seed)
    samples = int(T / every)
    out = np.zeros(samples)
    n, t, next_t, j = n0, 0.0, 0.0, 0
    while j < samples:
        kt = k + (k_end - k) * t / T
        u = n / N
        b = N * r * u
        d = N * (r * u * u / kt + u * u / (1 + u * u))
        total = b + d
        if total <= 0:
            out[j:] = 0.0
            break
        t += -np.log(np.random.random()) / total
        while j < samples and t >= next_t:
            out[j] = u
            j += 1
            next_t += every
        if np.random.random() < b / total:
            n += 1
        else:
            n -= 1
    return out


@njit
def budworm_switch(n0, target, N, r, k, reps, tmax, seed):
    np.random.seed(seed)
    out = np.zeros(reps)
    up = target > n0
    for i in range(reps):
        n, t = n0, 0.0
        while t < tmax and ((up and n < target) or ((not up) and n > target)):
            u = n / N
            b = N * r * u
            d = N * (r * u * u / k + u * u / (1 + u * u)) if n > 1 else 0.0   # reflecting at 1, as the exact chain
            t += -np.log(np.random.random()) / (b + d)
            n += 1 if np.random.random() < b / (b + d) else -1
        out[i] = t
    return out


def part_budworm(rng, quick):
    r, k = 0.55, 8.0
    f = lambda u, k=k: r * u * (1 - u / k) - u * u / (1 + u * u)
    fold = lambda u: 2 * u ** 3 / (1 + u * u) ** 2 - r
    u1, u2 = brentq(fold, 1.01, np.sqrt(3)), brentq(fold, np.sqrt(3), 6.0)
    k_up, k_down = 2 * u1 ** 3 / (u1 ** 2 - 1), 2 * u2 ** 3 / (u2 ** 2 - 1)
    low = brentq(f, 0.2, 1.2)
    mid = brentq(f, 1.2, 3.0)
    high = brentq(f, 3.0, 8.0)
    res = {"r": r, "k": k, "folds": [k_down, k_up], "states": [low, mid, high]}
    # hysteresis: the capacity swept slowly up and down, for several habitat sizes and two sweep times
    sweeps = []
    grid = np.linspace(3, 20, 341)
    for N in (100, 1000, 10000) + (() if quick else (100000,)):
        for T in (2000.0, 20000.0):
            n0 = int(N * brentq(lambda u: f(u, 3.0), 0.2, 3.0))
            up = budworm_events(n0, N, r, 3.0, T, T / 340, int(rng.integers(2 ** 31)), 20.0)
            n1 = int(N * brentq(lambda u: f(u, 20.0), 3.0, 30.0))
            dn = budworm_events(n1, N, r, 20.0, T, T / 340, int(rng.integers(2 ** 31)), 3.0)
            ks_up = 3 + 17 * np.arange(len(up)) / len(up)     # the capacity at each sample
            ks_dn = 20 - 17 * np.arange(len(dn)) / len(dn)
            thr = np.sqrt(low * high)
            jump_up = float(ks_up[np.argmax(up > thr)]) if (up > thr).any() else None
            jump_dn = float(ks_dn[np.argmax(dn < thr)]) if (dn < thr).any() else None
            sweeps.append({"N": N, "T": T, "jump_up": jump_up, "jump_down": jump_dn,
                           "k": ks_up[::2].tolist() if (N == 10000 and T == 20000.0) else None,
                           "up": up[::2].tolist() if (N == 10000 and T == 20000.0) else None,
                           "down": dn[::2].tolist() if (N == 10000 and T == 20000.0) else None})
    res["sweeps"] = sweeps
    # the time to leave a state in the bistable range: exact from the rates, simulated where it is short enough
    rows = []
    for N in (10, 15, 20, 30, 40, 60, 80, 100, 150, 200, 300):
        birth = lambda n, N=N: N * r * (n / N)
        death = lambda n, N=N: N * (r * (n / N) ** 2 / k + (n / N) ** 2 / (1 + (n / N) ** 2))
        n_low, n_mid, n_high = int(round(low * N)), int(round(mid * N)), int(round(high * N))
        if not (0 < n_low < n_mid < n_high):
            continue
        # chains with a reflecting lower end at n = 1 and a large upper end: extinction is excluded by the reflection
        up = exit_time(birth, death, n_low, n_mid)
        dn = exit_time(birth, death, n_high, n_mid, N_max=int(N * k * 1.6))
        row = {"N": N, "ln_tau_up": float(np.log(up)), "ln_tau_down": float(np.log(dn))}
        for name, s, tgt, tau in (("up", n_low, n_mid, up), ("down", n_high, n_mid, dn)):
            if tau < 3e3:
                t = budworm_switch(s, tgt, N, r, k, 200, 1e6, int(rng.integers(2 ** 31)))
                row[f"ln_tau_{name}_sim"] = float(np.log(t.mean()))
                row[f"ln_tau_{name}_sim_se"] = float(t.std() / t.mean() / np.sqrt(len(t)))
        rows.append(row)
    res["switching"] = rows
    trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    S_up = float(trap(np.log(r * np.linspace(low, mid, 20001) / (r * np.linspace(low, mid, 20001) ** 2 / k
                 + np.linspace(low, mid, 20001) ** 2 / (1 + np.linspace(low, mid, 20001) ** 2))), np.linspace(low, mid, 20001)))
    x = np.linspace(mid, high, 20001)
    S_dn = float(trap(np.log((r * x) / (r * x * x / k + x * x / (1 + x * x))), x))
    res["action_up"], res["action_down"] = -S_up, S_dn
    big = [row for row in rows if row["N"] >= 100]
    res["action_up_from_exact"] = float(np.polyfit([q["N"] for q in big], [q["ln_tau_up"] for q in big], 1)[0])
    res["action_down_from_exact"] = float(np.polyfit([q["N"] for q in big], [q["ln_tau_down"] for q in big], 1)[0])
    # constant noise, as a specification has it: barrier over the noise intensity at the state, D = (b + d) / (2N) = b / N
    xl, xh = np.linspace(low, mid, 20001), np.linspace(mid, high, 20001)
    V_up = float(-trap(f(xl), xl))
    V_dn = float(trap(f(xh), xh))
    res["constant_noise_action_up"] = V_up / (r * low)
    res["constant_noise_action_down"] = V_dn / (r * high)
    return res


# ------------------------------------------------------------------------------------------------ closed results
@njit
def budworm_jump(n0, stop, N, r, k0, k1, T, seed):
    """Capacity at which a sweep from k0 to k1 over time T first brings the population to the count stop."""
    np.random.seed(seed)
    up = stop > n0
    n, t = n0, 0.0
    while t < T and ((up and n < stop) or ((not up) and n > stop)):
        kt = k0 + (k1 - k0) * t / T
        u = n / N
        b = N * r * u
        d = N * (r * u * u / kt + u * u / (1 + u * u))
        t += -np.log(np.random.random()) / (b + d)
        n += 1 if np.random.random() < b / (b + d) else -1
    return k0 + (k1 - k0) * min(t, T) / T


def part_closed(rng, quick):
    res = {}
    # Kuramoto: dr/dt = (K/2 - Delta) r - (K/2) r^3 from all phases equal, fitted over the window used for the agents.
    # The cubic term makes the fitted rate exceed Delta - K/2 by about (K/2) r^2.
    delta, rel = 0.5, []
    for K in (0.2, 0.5, 0.7, 0.8, 0.9):
        t = np.arange(0.0, 60.0, 0.2)
        sol = solve_ivp(lambda t, y: [(K / 2 - delta) * y[0] - K / 2 * y[0] ** 3], (0, 60.0), [1.0], t_eval=t, rtol=1e-10,
                        atol=1e-13, method="LSODA")
        ok = (sol.y[0] < 0.25) & (sol.y[0] > 0.03)
        rel.append({"K": K, "rate": float(-np.polyfit(t[ok], np.log(sol.y[0][ok]), 1)[0])})
    res["kuramoto_relaxation_closed"] = rel

    r, k0 = 0.55, 8.0
    f = lambda u, k: r * u * (1 - u / k) - u * u / (1 + u * u)
    fold = lambda u: 2 * u ** 3 / (1 + u * u) ** 2 - r
    u1, u2 = brentq(fold, 1.01, np.sqrt(3)), brentq(fold, np.sqrt(3), 6.0)
    k_up, k_down = 2 * u1 ** 3 / (u1 ** 2 - 1), 2 * u2 ** 3 / (u2 ** 2 - 1)
    low0 = brentq(lambda u: f(u, k0), 0.2, 1.2)
    high0 = brentq(lambda u: f(u, k0), 3.0, 8.0)
    thr = float(np.sqrt(low0 * high0))
    res["budworm_threshold"] = thr

    def branches(k):
        return brentq(lambda u: f(u, k), 0.2, u1), brentq(lambda u: f(u, k), u1, u2), brentq(lambda u: f(u, k), u2, 40.0)

    # the sweep of the equation itself: the delay past each fold, as the jump is read for the agents
    det = []
    for T in (2000.0, 20000.0, 200000.0):
        row = {"T": T}
        for name, ka, kb, u0 in (("up", 3.0, 20.0, brentq(lambda u: f(u, 3.0), 0.2, 3.0)),
                                 ("down", 20.0, 3.0, brentq(lambda u: f(u, 20.0), 3.0, 30.0))):
            event = lambda t, y: y[0] - thr
            event.terminal = True
            sol = solve_ivp(lambda t, y: [f(y[0], ka + (kb - ka) * t / T)], (0, T), [u0], events=event, method="LSODA",
                            rtol=1e-9, atol=1e-11, max_step=T / 500)
            row[f"jump_{name}"] = float(ka + (kb - ka) * sol.t_events[0][0] / T)
        det.append(row)
    res["budworm_sweep_equation"] = det
    res["budworm_delay_exponent_up"] = float(np.polyfit(np.log([1 / d["T"] for d in det]),
                                                        np.log([d["jump_up"] - k_up for d in det]), 1)[0])
    res["budworm_delay_exponent_down"] = float(np.polyfit(np.log([1 / d["T"] for d in det]),
                                                          np.log([k_down - d["jump_down"] for d in det]), 1)[0])

    # the capacity at which the two states are equally occupied. With constant noise of any strength it is where the
    # potential difference vanishes; for birth and death it is where the integral of ln(b/d) vanishes.
    def constant_noise(k):
        a, _, c = branches(k)
        x = np.linspace(a, c, 20001)
        return float(TRAP(f(x, k), x))

    def births_deaths(k):
        a, _, c = branches(k)
        x = np.linspace(a, c, 20001)
        return float(TRAP(np.log(r / (r * x / k + x / (1 + x * x))), x))

    lo, hi = k_down + 1e-3, k_up - 1e-3
    res["budworm_equal_occupation"] = {"constant_noise": brentq(constant_noise, lo, hi), "birth_death": brentq(births_deaths, lo, hi)}

    # the same from the exact stationary distribution of the chain at finite N, with states divided at the middle state
    def log_weights(N, k):
        top = int(N * k * 1.6)
        n = np.arange(1, top + 1)
        u = n / N
        log_ratio = np.log(N * r * u[:-1]) - np.log(N * (r * u[1:] ** 2 / k + u[1:] ** 2 / (1 + u[1:] ** 2)))
        lp = np.concatenate(([0.0], np.cumsum(log_ratio)))
        cut = int(N * branches(k)[1])
        return logsumexp(lp[n <= cut]) - logsumexp(lp[n > cut])

    res["budworm_equal_occupation_chain"] = {str(N): brentq(lambda k: log_weights(N, k), lo, hi, xtol=1e-3) for N in (100, 300, 1000, 3000)}

    # the jump of a slow sweep at small N: activated. The chance to have left by k is 1 - exp(-int dk' (T/17) / tau(k')), with
    # tau the exact mean time to reach the threshold at fixed k. Repeated sweeps test the median it gives.
    def birth(N):
        return lambda n: N * r * (n / N)

    def death(N, k):
        return lambda n: N * (r * (n / N) ** 2 / k + (n / N) ** 2 / (1 + (n / N) ** 2))

    np.seterr(over="ignore")   # a time beyond the largest float is an infinite time here, with a rate of zero
    rows = []
    T = 20000.0
    for N, reps in (((100, 20),) if quick else ((100, 40), (300, 40), (1000, 24))):
        row = {"N": N, "T": T, "reps": reps}
        stop = int(np.ceil(thr * N))
        # upward
        ks = np.linspace(k_down + 0.02, k_up - 0.01, 120)
        rate = [1 / exit_time(birth(N), death(N, k), int(round(branches(k)[0] * N)), stop) for k in ks]
        H = np.concatenate(([0.0], np.cumsum(0.5 * (np.array(rate[1:]) + rate[:-1]) * np.diff(ks)))) * T / 17.0
        row["up_predicted"] = float(np.interp(np.log(2), H, ks)) if H[-1] > np.log(2) else None
        n0 = int(N * brentq(lambda u: f(u, 3.0), 0.2, 3.0))
        got = [budworm_jump(n0, stop, N, r, 3.0, 20.0, T, int(rng.integers(2 ** 31))) for _ in range(reps)]
        row["up_median"], row["up_sd"] = float(np.median(got)), float(np.std(got))
        # downward
        ks = np.linspace(k_up - 0.01, k_down + 0.02, 120)
        rate = [1 / exit_time(birth(N), death(N, k), int(round(branches(k)[2] * N)), stop, N_max=int(N * k * 1.6)) for k in ks]
        H = np.concatenate(([0.0], np.cumsum(0.5 * (np.array(rate[1:]) + rate[:-1]) * -np.diff(ks)))) * T / 17.0
        row["down_predicted"] = float(np.interp(np.log(2), H, ks)) if H[-1] > np.log(2) else None
        n1 = int(N * brentq(lambda u: f(u, 20.0), 3.0, 30.0))
        got = [budworm_jump(n1, stop, N, r, 20.0, 3.0, T, int(rng.integers(2 ** 31))) for _ in range(reps)]
        row["down_median"], row["down_sd"] = float(np.median(got)), float(np.std(got))
        rows.append(row)
    res["budworm_activated_jumps"] = rows
    return res


def part_fold_passage(rng, quick):
    """The jump of a sweep at large N, started on the branch at k = 8, at the rate of a full sweep from 3 to 20 in T."""
    r, k0 = 0.55, 8.0
    f = lambda u: r * u * (1 - u / k0) - u * u / (1 + u * u)
    low, high = brentq(f, 0.2, 1.2), brentq(f, 3.0, 8.0)
    thr = np.sqrt(low * high)
    rows = []
    for N, T, reps in ((10000, 2000.0, 8), (10000, 20000.0, 8), (100000, 2000.0, 5), (100000, 20000.0, 3)):
        stop = int(np.ceil(thr * N))
        up = [budworm_jump(int(low * N), stop, N, r, k0, 13.0 if T > 2000.0 else 14.0, T * (5.0 if T > 2000.0 else 6.0) / 17,
                           int(rng.integers(2 ** 31))) for _ in range(reps)]
        dn = [budworm_jump(int(high * N), stop, N, r, k0, 5.0, T * 3.0 / 17, int(rng.integers(2 ** 31))) for _ in range(reps)]
        rows.append({"N": N, "T": T, "reps": reps, "up": float(np.mean(up)), "up_sd": float(np.std(up)),
                     "down": float(np.mean(dn)), "down_sd": float(np.std(dn))})
    return {"jumps": rows}


def part_actions(rng, quick):
    """The exponent of the time to leave a state per unit N, in three descriptions of the same birth and death: the
    rates themselves, ln(b/d); a diffusion with the state-dependent intensity (b + d)/2N, 2(b - d)/(b + d); and a
    constant intensity, the value at the state being left."""
    def three(b, d, start, stop):
        x = np.linspace(start, stop, 40001)  # travel against the drift, so that each integral is positive
        return {"rates": float(TRAP(np.log(d(x) / b(x)), x)),
                "diffusion": float(TRAP(2 * (d(x) - b(x)) / (d(x) + b(x)), x)),
                "constant": float(TRAP(d(x) - b(x), x)) / float(b(np.array([start]))[0])}
    r, k = 0.55, 8.0
    f = lambda u: r * u * (1 - u / k) - u * u / (1 + u * u)
    low, mid, high = brentq(f, 0.2, 1.2), brentq(f, 1.2, 3.0), brentq(f, 3.0, 8.0)
    b = lambda u: r * u
    d = lambda u: r * u * u / k + u * u / (1 + u * u)
    res = {"budworm_up": three(b, d, low, mid), "budworm_down": three(b, d, high, mid)}
    f0, q = 2.0, 0.8
    x_star = (f0 * q - 1) / (f0 - 1)
    be = lambda x: f0 * q * x * (1 - x) / (1 + (f0 - 1) * x)
    de = lambda x: x * (f0 * x * (1 - q) + 1 - x) / (1 + (f0 - 1) * x)
    res["error"] = three(be, de, x_star, 1e-6)
    return res


TRAP = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
PARTS = {"kuramoto": part_kuramoto, "cycle": part_cycle, "error": part_error, "budworm": part_budworm, "closed": part_closed,
         "fold_passage": part_fold_passage, "actions": part_actions}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parts", nargs="*", default=list(PARTS), choices=list(PARTS))
    ap.add_argument("--quick", action="store_true", help="Smaller populations and shorter runs.")
    ap.add_argument("--seed", type=int, default=20261030)
    ap.add_argument("--out", default="results.json")
    args = ap.parse_args()
    out = Path(args.out)
    res = json.loads(out.read_text(encoding="utf-8")) if out.exists() else {}
    for part in args.parts:
        t0 = time.time()
        res[part] = PARTS[part](np.random.default_rng([args.seed, list(PARTS).index(part)]), args.quick)
        res[part]["seconds"] = round(time.time() - t0, 1)
        out.write_text(json.dumps(res, indent=1, default=float) + "\n", encoding="utf-8")
        print(part, "done in", res[part]["seconds"], "s", flush=True)


if __name__ == "__main__":
    main()
