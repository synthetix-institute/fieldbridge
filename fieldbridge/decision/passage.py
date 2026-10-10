"""Passage from a seed after a step, and the passage past the ghost of a fold.

A step (or a cue beyond the coercive field) makes a mode of a collective state unstable. The mode has d components:
d = 1 when one real eigenvalue leads (its symmetry is a reflection: a bend, a chiral excess, the difference of two
pools), d = 2 when a complex pair or a degenerate pair leads (a rotation: a transverse magnetization, a complex order
parameter). Along the mean trajectory s(t) = e^{Lambda(t)} (xi + ...), Lambda = int lambda dt, with xi a Gaussian seed
of d components and spread proportional to N^(-1/2). A threshold far from the start is crossed when
Lambda(tau) = ln(L/|xi|), so

    window:  Lambda(tau_90) - Lambda(tau_10) = ln q_d,  q_d = chi_d^{-1}(0.9) / chi_d^{-1}(0.1)
             (q_1 = 13.09, q_2 = 4.675);
    median:  Lambda(tau_50) rises by 1/2 per factor e in N;

for a constant rate, lambda Delta tau = ln q_d. A seed with mean z (one component, in units of its spread) selects the
sign with Phi(z) and narrows the window to the folded-normal ratio. The statistics need ln(L/sigma) >> 1; in practice
Lambda(tau_50) >= 2. These are the passage statistics of Haake, Haus and Glauber (1981); what a card certifies is the
rate, from the linearization of the model along its own mean trajectory.

Where a cue first removes a fold, the passage past its ghost takes pi / |a_h b (h - h_c)|^(1/2) from the normal form
s' = a_h (h - h_c) + b s^2.
"""
from __future__ import annotations

from math import log, pi, sqrt
from typing import Callable, Dict, Sequence

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi, norm

WINDOW = (0.1, 0.9)


def q_ratio(d: int, lo: float = 0.1, hi: float = 0.9) -> float:
    """chi_d^{-1}(hi) / chi_d^{-1}(lo): the ratio of seed magnitudes that bound the lo-hi window of passage times."""
    return float(chi.ppf(hi, d) / chi.ppf(lo, d))


def folded_ratio(z: float, lo: float = 0.1, hi: float = 0.9) -> float:
    """q_ratio for a one-component seed with mean z (folded normal F(x) = Phi(x - z) - Phi(-x - z))."""
    if z == 0:
        return q_ratio(1, lo, hi)
    F = lambda x, p: norm.cdf(x - z) - norm.cdf(-x - z) - p                      # noqa: E731
    return float(brentq(F, 0.0, abs(z) + 20, args=(hi,)) / brentq(F, 0.0, abs(z) + 20, args=(lo,)))


def window_law(d: int, z: float = 0.0) -> float:
    """ln q_d (or ln of the folded ratio for a biased one-component seed)."""
    return log(folded_ratio(z)) if (d == 1 and z) else log(q_ratio(d))


class Exponent:
    """Lambda(t) = int_0^t lambda dt on a grid, from the rate lambda(t) along the mean trajectory."""

    def __init__(self, t: np.ndarray, lam: np.ndarray):
        self.t = np.asarray(t, float)
        self.lam = np.asarray(lam, float)
        self.L = np.concatenate([[0.0], np.cumsum(0.5 * (self.lam[1:] + self.lam[:-1]) * np.diff(self.t))])

    def __call__(self, tau):
        return np.interp(tau, self.t, self.L)

    def rate(self, tau):
        return np.interp(tau, self.t, self.lam)


def quantiles(times: np.ndarray, boot: int = 400, seed: int = 0) -> Dict[str, float]:
    """10, 50 and 90% quantiles of the finite passage times, with bootstrap standard errors."""
    t = np.asarray(times, float)
    t = t[np.isfinite(t)]
    if t.size < 20:
        raise ValueError("fewer than 20 replicas crossed the threshold")
    ps = (WINDOW[0], 0.5, WINDOW[1])
    q = np.quantile(t, ps)
    rng = np.random.default_rng(seed)
    bs = np.array([np.quantile(rng.choice(t, t.size), ps) for _ in range(boot)])
    return {"q10": float(q[0]), "q50": float(q[1]), "q90": float(q[2]), "q10_se": float(bs[:, 0].std()),
            "q50_se": float(bs[:, 1].std()), "q90_se": float(bs[:, 2].std()), "window": float(q[2] - q[0]),
            "finite": int(t.size)}


def check(times: np.ndarray, ex: Exponent, d: int, z: float = 0.0, tol: float = 0.03) -> Dict[str, float]:
    """The measured window in Lambda against the law, with its bootstrap error (gate |dL - law| < 4 se + tol law)."""
    q = quantiles(times)
    dL = float(ex(q["q90"]) - ex(q["q10"]))
    se = float(np.hypot(ex.rate(q["q90"]) * q["q90_se"], ex.rate(q["q10"]) * q["q10_se"]))
    law = window_law(d, z)
    return {**q, "d": d, "z": z, "Lambda_window": dL, "Lambda_window_se": se, "law": law, "ratio": dL / law,
            "agrees": bool(abs(dL - law) < 4 * se + tol * law), "Lambda_median": float(ex(q["q50"])),
            "Lambda_median_se": float(ex.rate(q["q50"]) * q["q50_se"])}


def median_slope(sizes: Sequence[float], medians: Sequence[float], errors: Sequence[float]) -> Dict[str, float]:
    """Slope of Lambda(tau_50) against ln N by weighted least squares (law: 1/2)."""
    x = np.log(np.asarray(sizes, float))
    y = np.asarray(medians, float)
    W = 1 / np.maximum(np.asarray(errors, float), 1e-9) ** 2
    A = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(A.T @ (W[:, None] * A))
    beta = cov @ A.T @ (W * y)
    return {"slope": float(beta[1]), "slope_se": float(sqrt(cov[1, 1])), "law": 0.5}


def leading_dimension(J: np.ndarray, tol: float = 1e-6) -> int:
    """d: the number of eigenvalues that share the largest real part (a complex pair counts 2)."""
    ev = np.linalg.eigvals(J)
    top = ev.real.max()
    return int(np.sum(np.abs(ev.real - top) <= tol * max(1.0, abs(top))))


def fold_normal_form(F: Callable, jac: Callable, x_c: np.ndarray, h_c: float, dh: float = 1e-6,
                     dx: float = 1e-4) -> Dict[str, float]:
    """At a saddle-node (x_c, h_c) of x' = F(x, h): a_h = w . dF/dh and b = (1/2) w . D2F(v, v), with v, w the right and
    left null vectors (w . v = 1); the ghost passage for h beyond h_c takes pi / |a_h b (h - h_c)|^(1/2)."""
    J = jac(x_c, h_c)
    ev, V = np.linalg.eig(J)
    v = np.real(V[:, np.argmin(np.abs(ev))])
    evl, Wl = np.linalg.eig(J.T)
    w = np.real(Wl[:, np.argmin(np.abs(evl))])
    w = w / (w @ v)
    a_h = float(w @ (F(x_c, h_c + dh) - F(x_c, h_c - dh)) / (2 * dh))
    b = float(w @ (F(x_c + dx * v, h_c) + F(x_c - dx * v, h_c) - 2 * F(x_c, h_c)) / (2 * dx * dx))
    return {"a_h": a_h, "b": b, "coefficient": pi / sqrt(abs(a_h * b)), "eigenvalues": ev.real.tolist()}
