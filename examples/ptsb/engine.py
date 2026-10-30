"""Batched trajectories for many surfaces and parameters at once, on NumPy or (if available) torch on a GPU.

Both model surfaces of ptsb.py are the polynomial

    V = a0 x^3 + a1 x^2 + a2 y^2 + a3 x y^2 + a4 x y^4 + a5 y + a6 x y^3 + a7 x^2 y^2 + a8 x^2 y + a9 x y
        + a10 y^4,

so a batch carries eleven coefficients per trajectory together with its friction, temperature and the geometry of its
product basins (commit level, x_commit, the TS2 ridge line, the reactant boundary). The integrator and the
assignment rules are those of ptsb.run.
"""
from __future__ import annotations

from typing import Dict, List

import numpy as np

from ptsb import Collins, Minimal, Surface


def coeffs(s: Surface) -> np.ndarray:
    if isinstance(s, Minimal):
        return np.array([8.0 / 3.0, -4.0, 0.5 + 2.0 * s.d, -2.0, 1.0, s.eps, s.c, 0.0, 0.0, 0.0, -s.d])
    if isinstance(s, Collins):
        c1, c2, c3, c4, c5, c6 = s.k
        return np.array([s.c0 / 3.0, -s.c0 * s.alpha / 2.0, s.om2 / 2.0, -s.om2 * s.beta / 2.0 + c4, c1, 0.0, c6,
                         c2, c3, c5, 0.0])
    raise TypeError(s)


def V_poly(a, x, y):
    y2 = y * y
    return (a[0] * x ** 3 + a[1] * x * x + a[2] * y2 + a[3] * x * y2 + a[4] * x * y2 * y2 + a[5] * y
            + a[6] * x * y2 * y + a[7] * x * x * y2 + a[8] * x * x * y + a[9] * x * y + a[10] * y2 * y2)


def grad_poly(a, x, y):
    y2 = y * y
    Vx = (3 * a[0] * x * x + 2 * a[1] * x + a[3] * y2 + a[4] * y2 * y2 + a[6] * y2 * y + 2 * a[7] * x * y2
          + 2 * a[8] * x * y + a[9] * y)
    Vy = (2 * a[2] * y + 2 * a[3] * x * y + 4 * a[4] * x * y2 * y + a[5] + 3 * a[6] * x * y2 + 2 * a[7] * x * x * y
          + a[8] * x * x + a[9] * x + 4 * a[10] * y2 * y)
    return Vx, Vy


def geometry(s: Surface) -> np.ndarray:
    qr, nr = s.ridge_normal()
    return np.array([s.commit_level(), s.x_commit(), qr[0], qr[1], nr[0], nr[1], s.points()["TS1"][0] - 0.4])


def run_batch(points: List[Dict], dt: float = 0.005, check_every: int = 4, backend: str = "numpy", seed: int = 0,
              device: str = "cuda") -> List[np.ndarray]:
    """points: dicts with keys surface, init (x, y, px, py arrays), gamma, kT, t_max, rule ("first" or "final").
    Returns the outcome array of each point (+1/-1 product side, 0 reactant, 2 undecided), and the first-entry
    assignment as a second array: [(outcome, first), ...]. A trajectory leaving |x|, |y| < 5 (the open region behind
    TS1 where both surfaces fall without bound) counts as returned to the reactant side."""
    sizes = [len(p["init"][0]) for p in points]
    cat = lambda f: np.concatenate([np.full(n, f(p), float) if np.isscalar(f(p)) else np.asarray(f(p), float)
                                    for n, p in zip(sizes, points)])
    A = np.concatenate([np.repeat(coeffs(p["surface"])[:, None], n, axis=1) for n, p in zip(sizes, points)], axis=1)
    Gm = np.concatenate([np.repeat(geometry(p["surface"])[:, None], n, axis=1) for n, p in zip(sizes, points)], axis=1)
    x, y, px, py = (np.concatenate([np.asarray(p["init"][i], float) for p in points]) for i in range(4))
    gamma, kT, tmax = cat(lambda p: p["gamma"]), cat(lambda p: p["kT"]), cat(lambda p: p["t_max"])
    final = cat(lambda p: 1.0 if p.get("rule", "first") == "final" else 0.0) > 0
    c1 = np.exp(-gamma * dt)
    c2 = np.sqrt((1.0 - c1 * c1) * kT)
    steps = int(round(tmax.max() / dt))
    if backend == "torch":
        import torch
        dev = torch.device(device)
        T = lambda arr: torch.as_tensor(arr, dtype=torch.float64, device=dev)
        A, Gm, x, y, px, py, c1, c2, tmax = (T(v) for v in (A, Gm, x, y, px, py, c1, c2, tmax))
        final = torch.as_tensor(final, device=dev)
        gen = torch.Generator(device=dev).manual_seed(seed)
        normal = lambda: torch.randn(x.shape, generator=gen, dtype=torch.float64, device=dev)
        full = lambda v: torch.full(x.shape, v, dtype=torch.int64, device=dev)
        where, sign, anyf = torch.where, torch.sign, lambda m: bool(m.any())
        to_np = lambda v: v.cpu().numpy()
    else:
        rng = np.random.default_rng(seed)
        normal = lambda: rng.standard_normal(x.shape)
        full = lambda v: np.full(x.shape, v, np.int64)
        where, sign, anyf = np.where, np.sign, lambda m: bool(m.any())
        to_np = lambda v: v
    noisy = anyf(c2 > 0)
    level, xc, qx, qy, nx, ny, xback = (Gm[i] for i in range(7))
    out, first = full(2), full(2)
    alive = full(1) > 0
    gx, gy = grad_poly(A, x, y)
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
        gx, gy = grad_poly(A, x, y)
        px = px - 0.5 * dt * gx
        py = py - 0.5 * dt * gy
        if (k + 1) % check_every == 0:
            t = (k + 1) * dt
            inside = (x > xc) & (V_poly(A, x, y) < level)
            sd = sign((x - qx) * nx + (y - qy) * ny)
            sd = sd.long() if backend == "torch" else sd.astype(np.int64)
            newly = alive & inside & (first == 2)
            first = where(newly, sd, first)
            back = alive & ((x < xback) | (abs(x) > 5.0) | (abs(y) > 5.0))
            first = where(back & (first == 2), full(0), first)
            out = where(back, full(0), out)
            done_first = alive & inside & ~final
            out = where(done_first, sd, out)
            expire = alive & ~back & ~done_first & (t >= tmax - 1e-9)
            out = where(expire & final & inside, sd, out)
            alive = alive & ~back & ~done_first & ~expire
            if not anyf(alive):
                break
            # freeze finished trajectories where they are
            px = where(alive, px, px * 0)
            py = where(alive, py, py * 0)
            gx = where(alive, gx, gx * 0)
            gy = where(alive, gy, gy * 0)
    out, first = to_np(out), to_np(first)
    res, i = [], 0
    for n in sizes:
        res.append((out[i:i + n], first[i:i + n]))
        i += n
    return res
