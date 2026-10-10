"""The onset of synchronization from the loss of stability of incoherence.

Phase oscillators with natural frequencies of density g(omega), coupled all to all through the first harmonic of
H(phi) = a1 cos phi + b1 sin phi,

    dtheta_i/dt = omega_i + (K/N) sum_j H(theta_j - theta_i),

have an incoherent state whose perturbations grow as e^{lambda t} with

    1 = (K/2) (b1 - i a1) int g(omega) / (lambda + i omega) d omega.

A mode turns unstable where lambda = epsilon - i Omega with epsilon -> 0+:

    2/K = (b1 - i a1) G(Omega),   G(Omega) = pi g(Omega) - i PV int g(omega)/(omega - Omega) d omega,

two real conditions for K_c and the frequency Omega of the emerging rhythm. For a1 = 0 and a symmetric unimodal g,
K_c = 2/(pi g(0) b1) (Kuramoto 1975); a1 != 0 moves K_c and Omega (Sakaguchi and Kuramoto 1986); b1 <= 0 has no onset.
Above K_c the coherent perturbation grows at the rate mu, the real part of the unstable root.
"""
from __future__ import annotations

from math import pi
from typing import Callable, Dict, Optional

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, fsolve


def G(g: Callable[[float], float], Omega: float, support=(-np.inf, np.inf), scale: float = 1.0) -> complex:
    """pi g(Omega) - i PV int g(w)/(w - Omega) dw (Cauchy weight on a window around Omega, plain quadrature outside)."""
    lo, hi = support
    width = 50.0 * scale
    a, b = max(lo, Omega - width), min(hi, Omega + width)
    pv = 0.0
    if a < Omega < b:
        pv += quad(g, a, b, weight="cauchy", wvar=Omega, limit=400)[0]
    else:
        pv += quad(lambda w: g(w) / (w - Omega), a, b, limit=400)[0]
    if lo < a:
        pv += quad(lambda w: g(w) / (w - Omega), lo, a, limit=400)[0]
    if b < hi:
        pv += quad(lambda w: g(w) / (w - Omega), b, hi, limit=400)[0]
    return complex(pi * g(Omega), -pv)


def critical_coupling(g: Callable[[float], float], a1: float, b1: float, centre: float = 0.0, scale: float = 1.0,
                      support=(-np.inf, np.inf), span: float = 6.0, grid: int = 241) -> Optional[Dict[str, float]]:
    """K_c and Omega: the roots of Im[(b1 - i a1) G(Omega)] = 0 in centre +- span*scale with positive real part, the
    smallest K among them. None when no root gives a positive K (no onset at any coupling). A root must lie where the
    density is positive: outside the band of frequencies g = 0 and Im[...] = -b1 PV, which for b1 = 0 vanishes at the
    edges of the band without a mode that grows (the coupling then only shifts the frequencies)."""
    c = complex(b1, -a1)
    g_ref = g(centre)
    f = lambda O: (c * G(g, O, support, scale)).imag             # noqa: E731
    Os = np.linspace(centre - span * scale, centre + span * scale, grid)
    vals = np.array([f(O) for O in Os])
    best: Optional[Dict[str, float]] = None
    for k in range(len(Os) - 1):
        if vals[k] == 0 or vals[k] * vals[k + 1] < 0:
            O = Os[k] if vals[k] == 0 else brentq(f, Os[k], Os[k + 1], xtol=1e-12 * scale)
            if g(O) <= 1e-6 * g_ref:
                continue
            re = (c * G(g, O, support, scale)).real
            if re > 0:
                K = 2.0 / re
                if best is None or K < best["K_c"]:
                    best = {"K_c": K, "Omega": O}
    return best


def growth_rate(g: Callable[[float], float], a1: float, b1: float, K: float, centre: float = 0.0, scale: float = 1.0,
                support=(-np.inf, np.inf), guess: Optional[Dict[str, float]] = None) -> Dict[str, float]:
    """The unstable root lambda = mu - i Omega (mu > 0) at the coupling K."""
    c = complex(b1, -a1)
    lo, hi = support
    lo, hi = max(lo, centre - 60 * scale), min(hi, centre + 60 * scale)

    def F(v):
        mu, Om = v
        mu = abs(mu)
        re = quad(lambda w: g(w) * mu / (mu ** 2 + (w - Om) ** 2), lo, hi, points=[Om], limit=400)[0]
        im = quad(lambda w: -g(w) * (w - Om) / (mu ** 2 + (w - Om) ** 2), lo, hi, points=[Om], limit=400)[0]
        z = 0.5 * K * c * complex(re, im) - 1.0
        return [z.real, z.imag]

    x0 = [guess["mu"], guess["Omega"]] if guess else [0.5 * K * b1, centre + 0.5 * K * a1]
    sol, info, ok, msg = fsolve(F, x0, full_output=True, xtol=1e-12)
    res = F(sol)
    return {"mu": float(abs(sol[0])), "Omega": float(sol[1]), "residual": float(np.hypot(*res)), "ok": bool(ok == 1)}


def lorentzian(gamma: float, centre: float = 0.0):
    return lambda w: gamma / pi / ((w - centre) ** 2 + gamma ** 2)


def gaussian(sigma: float, centre: float = 0.0):
    return lambda w: float(np.exp(-0.5 * ((w - centre) / sigma) ** 2) / (sigma * np.sqrt(2 * pi)))


def uniform(width: float, centre: float = 0.0):
    return lambda w: (1.0 / (2 * width)) if abs(w - centre) < width else 0.0
