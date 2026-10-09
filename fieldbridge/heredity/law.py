"""Inheritance through a threshold: the law.

A daughter starts below the threshold of a pitchfork with the order phi0 of its parent and is carried through the
threshold by growth. In the normal form

    dphi/dt = a mu(t) phi - b phi^3 + sqrt(2D) xi,    mu(t) = -mu0 + r t,

it keeps the parent's sign with probability

    P = Phi(phi_c / sigma_c),   sigma_c^2 = 2D sqrt(pi/(a r)) Phi(a mu0 sqrt(2/(a r))),

where phi_c is the deterministic order at the crossing time mu0/r and sigma_c the noise accumulated by the critical mode
in the linear passage. With a sustained bias h in place of an initial order the same calculation gives the swept-write
law P = Phi(pi^(1/4) h / (D^(1/2) (a r)^(1/4))) of the memory module.

phi_c must be the order that the body's own equations carry to the crossing: the linear decay phi0 e^(-a mu0^2/(2r)) of
the normal form misses the cubic term (off by up to 0.21 over 120 conditions of the normal form) and, in a body whose
critical eigenvalue is not linear in the size through a deep dip, the local normal form misses the decay itself.
"""
from __future__ import annotations

from math import erf, exp, pi, sqrt

import numpy as np
from scipy.integrate import solve_ivp


def Phi(x: float) -> float:
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))


def sigma_c(D: float, mu0: float, r: float, a: float = 1.0) -> float:
    """Standard deviation of the noise accumulated by the critical mode from the start to the end of the passage, in
    units of phi at the crossing time: 2D int_{-t0}^{inf} exp(-a r s^2) ds with t0 = mu0/r."""
    ar = a * r
    return sqrt(2.0 * D * sqrt(pi / ar) * Phi(a * mu0 * sqrt(2.0 / ar)))


def phi_at_crossing(phi0: float, mu0: float, r: float, a: float = 1.0, b: float = 1.0) -> float:
    """phi at the crossing time mu0/r of the deterministic normal form (with the cubic term)."""
    tc = mu0 / r
    if tc <= 0:
        return phi0
    sol = solve_ivp(lambda t, y: [a * (-mu0 + r * t) * y[0] - b * y[0] ** 3], (0.0, tc), [phi0],
                    rtol=1e-11, atol=1e-15)
    return float(sol.y[0, -1])


def probability(phi0: float, mu0: float, r: float, D: float, a: float = 1.0, b: float = 1.0,
                cubic: bool = True) -> float:
    """P of keeping the parent's sign (phi0 > 0) through the passage, in the normal form."""
    if D <= 0:
        return 1.0 if phi0 > 0 else 0.5
    phic = phi_at_crossing(phi0, mu0, r, a, b) if cubic else phi0 * exp(-a * mu0 ** 2 / (2 * r))
    return Phi(phic / sigma_c(D, mu0, r, a))


def probabilities(sign: np.ndarray, phic: np.ndarray, sigma: float) -> np.ndarray:
    """Phi(sign phi_c / sigma) for each daughter; without noise the sign at the crossing decides (1/2 at phi_c = 0)."""
    x = np.asarray(sign) * np.asarray(phic)
    if sigma <= 0:
        return np.where(x > 0, 1.0, np.where(x < 0, 0.0, 0.5))
    return np.array([Phi(v) for v in x / sigma])


def simulate_normal_form(phi0: float, mu0: float, r: float, D: float, n: int, rng, a: float = 1.0, b: float = 1.0,
                         dt: float = 2e-3, beyond: float = 1.0) -> float:
    """Fraction of n daughters that keep the parent's sign: Euler-Maruyama from t = 0 until mu = beyond, well past the
    threshold, where the order has regrown and its sign no longer changes."""
    phi = np.full(n, float(phi0))
    t, T = 0.0, (mu0 + beyond) / r
    amp = sqrt(2.0 * D * dt)
    while t < T:
        phi += (a * (-mu0 + r * t) * phi - b * phi ** 3) * dt + amp * rng.standard_normal(n)
        t += dt
    return float(np.mean(phi > 0))
