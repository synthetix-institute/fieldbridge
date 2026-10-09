"""Lineages through growth and division, and the law of inheritance through a threshold (letters K and L).

A body whose order exists only above the size L_c grows from L_div/2 to L_div and divides; one daughter is followed.
The daughter starts below L_c, its order decays, and growth carries it back through the pitchfork. It inherits when the
sign of its order at its own division equals its parent's. The order is w.(q - q_sym), with w the left critical vector.

The law in three versions:
    normal form   law.probability with a linear ramp dL/dt at L_c and the cubic decay of the reduced equation;
    body          the order at the crossing from the body's noise-free equations along the actual growth, with the
                  noise of the linear passage, law.sigma_c;
    history       the same order with the noise accumulated along the critical eigenvalue's own history,
                  sigma^2 = int 2 D_s(L(t)) exp(-2 (Lambda(t) - Lambda(t_c))) dt.
The body version is the law of the module: the normal form fails where the dip is deep (the critical eigenvalue is not
linear in the size there). Each row reports the conditions of the law: Lambda = a ramp / (b D_s) >> 1 (small noise in
the crossing window) and ln G = int lambda dt over one generation > 0 (the order survives a generation without noise).
"""
from __future__ import annotations

from math import sqrt
from typing import Dict, Optional

import numpy as np

from . import law


def noise_s(lin, red: Dict, L: float, noise: float = 1.0) -> float:
    """D_s = w^T D w at the symmetric state of size L."""
    q = lin.body.steady(red["q_sym"], L)
    D = np.asarray(lin.body.diffusion(q, L), float)
    w = red["w"]
    return float(noise * (w @ D @ w if D.ndim == 2 else np.sum(w ** 2 * D)))


def order(red: Dict, Q: np.ndarray) -> np.ndarray:
    return (np.atleast_2d(Q) - red["q_sym"]) @ red["w"]


def history(lin, red: Dict, L0: float, L1: float, points: int = 120):
    """Along the growth from L0 to L1: times, sizes, the critical eigenvalue lambda at the symmetric state, its integral
    Lambda(t), and D_s = w^T D w (noise of the specification). Computed once per interval and kept in red."""
    cache = red.setdefault("_histories", {})
    key = (round(L0, 12), round(L1, 12), points)
    if key in cache:
        return cache[key]
    ts = np.linspace(0.0, lin.growth.duration(L0, L1), points)
    Ls = lin.growth.size(L0, ts)
    w, v = red["w"], red["v"]
    q = red["q_sym"]
    lam, Ds = np.empty(points), np.empty(points)
    for k, L in enumerate(Ls):
        q = lin.body.steady(q, float(L))
        ev, V = np.linalg.eig(lin.body.jac(q, float(L)))
        Vr = np.real(V)
        lam[k] = ev[int(np.argmax(np.abs(Vr.T @ v) / np.linalg.norm(Vr, axis=0).clip(1e-300)))].real
        D = np.asarray(lin.body.diffusion(q, float(L)), float)
        Ds[k] = w @ D @ w if D.ndim == 2 else np.sum(w ** 2 * D)
    Lam = np.concatenate([[0.0], np.cumsum(0.5 * (lam[1:] + lam[:-1]) * np.diff(ts))])
    cache[key] = (ts, Ls, lam, Lam, Ds)
    return cache[key]


def eigenvalue_history(lin, red: Dict, L0: float, L1: float, points: int = 120):
    ts, Ls, lam, Lam, _ = history(lin, red, L0, L1, points)
    return ts, Ls, lam, Lam


def ln_gain(lin, red: Dict, L_div: float) -> float:
    """int lambda dt over one generation, L_div/2 -> L_div: the order survives without noise when it is positive."""
    return float(history(lin, red, L_div / 2, L_div)[3][-1])


def sigma_history(lin, red: Dict, L0: float, L_end: float, noise: float = 1.0) -> float:
    ts, Ls, lam, Lam, Ds = history(lin, red, L0, L_end)
    Lam_c = np.interp(lin.growth.duration(L0, red["L_c"]), ts, Lam)
    integrand = 2 * noise * Ds * np.exp(-2 * (Lam - Lam_c))
    return sqrt(float(np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(ts))))


def noise_free_lineage(lin, red: Dict, L_div: float, generations: int = 3, points: int = 60,
                       start_order: Optional[float] = None, seed: int = 0):
    """One lineage without noise. The first parent starts at the threshold with the order start_order (default: the
    noise amplitude sigma_c of the linear passage, as a noisy parent leaves the crossing) and grows to L_div; then
    `generations` daughters follow. Returns the times, sizes and orders at `points` instants per generation and the
    state of each daughter at birth."""
    rng = np.random.default_rng(seed)
    Lc = red["L_c"]
    if start_order is None:
        Ds = noise_s(lin, red, Lc)
        start_order = law.sigma_c(Ds, Lc - L_div / 2, lin.growth.ramp(Lc), red["a"]) if Ds > 0 else 0.05
        start_order = max(start_order, 1e-6)
    q = (lin.body.steady(red["q_sym"], Lc) + start_order * red["v"])[None, :]
    ts, Ls, ss, births = [], [], [], []
    t0 = 0.0
    for g in range(generations + 1):
        L0 = Lc if g == 0 and L_div / 2 < Lc else L_div / 2
        dur = lin.growth.duration(L0, L_div)
        grid = lin.growth.size(L0, np.linspace(0.0, dur, points + 1))
        for k in range(points):
            ts.append(t0 + dur * k / points)
            Ls.append(float(grid[k]))
            ss.append(float(order(red, q)[0]))
            q = lin.body.grow(q, float(grid[k]), float(grid[k + 1]), rng, noise=0.0)
        t0 += dur
        ts.append(t0)
        Ls.append(L_div)
        ss.append(float(order(red, q)[0]))
        q = lin.body.divide(q, rng, noise=0.0)
        births.append(q[0].copy())
    return {"t": ts, "L": Ls, "order": ss, "births": births, "start_order": start_order}


def law_for_daughter(lin, red: Dict, q_born: np.ndarray, L_div: float, noise: float = 1.0) -> Dict:
    """The law for a daughter born in the state q_born: the order at birth and at the crossing, sigma (linear passage
    and history), and P in its three versions."""
    rng = np.random.default_rng(0)
    L0, Lc = L_div / 2, red["L_c"]
    s0 = float(order(red, q_born)[0])
    phic = float(order(red, lin.body.grow(np.atleast_2d(q_born), L0, Lc, rng, noise=0.0))[0])
    ramp, mu0 = lin.growth.ramp(Lc), Lc - L0
    Ds = noise_s(lin, red, Lc, noise)
    sig_lin = law.sigma_c(Ds, mu0, ramp, red["a"]) if Ds > 0 else 0.0
    sig_hist = sigma_history(lin, red, L0, L_div, noise) if Ds > 0 else 0.0
    sign = np.sign(s0)
    return {"order_at_birth": abs(s0), "order_at_crossing": float(sign * phic), "D_s": Ds,
            "sigma_linear": sig_lin, "sigma_history": sig_hist,
            "P_body": float(law.probabilities([sign], [phic], sig_lin)[0]),
            "P_history": float(law.probabilities([sign], [phic], sig_hist)[0]),
            "P_normal_form": law.probability(abs(s0), mu0, ramp, Ds, a=red["a"], b=red["b"]) if s0 != 0 else 0.5,
            "Lambda": red["a"] * ramp / (red["b"] * Ds) if red["b"] > 0 and Ds > 0 else float("inf")}


def condition(lin, red: Dict, L_div: float, noise: float = 1.0, n: int = 2000, generations: int = 4, burn: int = 1,
              seed: int = 1) -> Dict:
    """Lineages at one condition and the law averaged over the daughters' states at birth."""
    rng = np.random.default_rng(seed)
    body = lin.body
    Q = np.repeat(red["q_sym"][None, :], n, 0) + lin.start_spread * rng.standard_normal((n, body.n))
    L0 = L_div / 2
    signs, kids = [], []
    for _ in range(generations):
        Q = body.grow(Q, L0, L_div, rng, noise=noise)
        signs.append(np.sign(order(red, Q)))
        Q = body.divide(Q, rng, noise=noise)
        kids.append(Q.copy())
    sign_p = np.concatenate(signs[burn:-1])
    sign_d = np.concatenate(signs[burn + 1:])
    born = np.concatenate(kids[burn:-1])
    measured = float(np.mean(sign_d == sign_p))
    s0 = order(red, born)
    Lc, a, b = red["L_c"], red["a"], red["b"]
    mu0, ramp = Lc - L0, lin.growth.ramp(Lc)
    Ds = noise_s(lin, red, Lc, noise)
    cache: Dict[float, float] = {}
    nf = []
    for x in np.abs(s0):
        key = round(float(x), 9)
        if key not in cache:
            cache[key] = law.probability(key, mu0, ramp, Ds, a=a, b=b) if key > 0 else 0.5
        nf.append(cache[key])
    nf = np.where(np.sign(s0) == sign_p, nf, 1 - np.array(nf))
    phic = order(red, body.grow(born, L0, Lc, rng, noise=0.0))
    sig_lin = law.sigma_c(Ds, mu0, ramp, a) if Ds > 0 else 0.0
    sig_hist = sigma_history(lin, red, L0, L_div, noise) if Ds > 0 else 0.0
    stderr = sqrt(max(measured * (1 - measured), 1e-12) / len(sign_p))
    return {"L_div": L_div, "noise": noise, "lineages": n, "transfers": int(len(sign_p)),
            "measured": measured, "stderr": stderr,
            "law_normal_form": float(np.mean(nf)),
            "law_body": float(np.mean(law.probabilities(sign_p, phic, sig_lin))),
            "law_history": float(np.mean(law.probabilities(sign_p, phic, sig_hist))),
            "sigma_linear": sig_lin, "sigma_history": sig_hist,
            "mean_abs_order_at_birth": float(np.mean(np.abs(s0))),
            "mean_abs_order_at_crossing": float(np.mean(np.abs(phic))),
            "Lambda": a * ramp / (b * Ds) if b > 0 and Ds > 0 else float("inf"),
            "lnG_generation": ln_gain(lin, red, L_div)}


def agrees(row: Dict, key: str = "law_body") -> bool:
    """The gate of the module: |P - P_law| < 4 stderr + 0.02."""
    return abs(row["measured"] - row[key]) < 4 * row["stderr"] + 0.02
