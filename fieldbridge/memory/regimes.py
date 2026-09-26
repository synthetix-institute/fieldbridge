"""Information about a write in the three regimes of kappa, and why it is not lost at kappa < 0.

For the linear stage near a state, dx = -kappa x dt + sqrt(2D) dW, a write that places the state at distance s0 (a
Gaussian write of variance s0^2) leaves the information

    I(t) = (1/2) ln[1 + (kappa s0^2 / D) / (exp(2 kappa t) - 1)]

about itself at time t. One expression covers the three regimes:
  kappa > 0   exponential loss (Law 1): the contraction dilutes the write in the stationary noise;
  kappa -> 0  I -> (1/2) ln(1 + s0^2 / (2 D t)), algebraic loss (Law 2);
  kappa < 0   I -> (1/2) ln(1 + W), W = |kappa| s0^2 / D: the loss stops. The information still only decreases (the
              evolution cannot create it), but the expansion outruns the noise: only the noise injected during the
              first ~1/|kappa|, near the ridge, competes with the write, and what survives it is frozen in. W is the
              number that sets the write law.

The binary version (write +-s0, measure the sign) is simulated for the three cases, the last with the saturating
pitchfork, where the wells then retain what was decided near the ridge. There the sign is set during the linear
stage: x e^{-|kappa| t} tends to s0 plus a Gaussian of variance D/|kappa|, so the sign is kept with probability
Phi(W^{1/2}), and the binary information stops at ln 2 - H(Phi(W^{1/2})).
"""
from __future__ import annotations

from typing import Dict, Sequence

import numpy as np


def gaussian_information(kappa: float, s0: float, D: float, t: np.ndarray) -> np.ndarray:
    t = np.asarray(t, float)
    if abs(kappa) < 1e-12:
        snr = s0 ** 2 / (2 * D * t)
    else:
        snr = (kappa * s0 ** 2 / D) / np.expm1(2 * kappa * t)
    return 0.5 * np.log1p(snr)


def _binary_information(p: np.ndarray) -> np.ndarray:
    """ln 2 - H(p) for a symmetric binary channel with accuracy p (nats)."""
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return np.log(2) + p * np.log(p) + (1 - p) * np.log(1 - p)


def simulate(kappa: float, s0: float, D: float, times: Sequence[float], rng, n: int = 20000, dt: float = 1e-3,
             saturate: bool = True) -> Dict[str, object]:
    """Binary write +-s0, observable sign(x). kappa < 0 uses the pitchfork dx = (|kappa| x - x^3) dt (saturating) if
    `saturate`, the linear instability otherwise."""
    times = np.asarray(times, float)
    sign = rng.choice([-1.0, 1.0], size=n)
    x = sign * s0
    amp = np.sqrt(2 * D * dt)
    out, t, k = [], 0.0, 0
    for target in times:
        while t < target - 1e-12:
            drift = -kappa * x - (x ** 3 if (kappa < 0 and saturate) else 0.0)
            x = x + dt * drift + amp * rng.standard_normal(n)
            t += dt
        out.append(float(np.mean(np.sign(x) == sign)))
    p = np.array(out)
    return {"times": times.tolist(), "accuracy": p.tolist(), "information": _binary_information(p).tolist()}


def binary_information_linear(kappa: float, s0: float, D: float, t: np.ndarray) -> np.ndarray:
    """Exact accuracy of the sign of the observable for the linear stage, turned into information."""
    from math import erf, sqrt
    t = np.asarray(t, float)
    if abs(kappa) < 1e-12:
        var = 2 * D * t
        mean = s0 * np.ones_like(t)
    else:
        var = (D / kappa) * (-np.expm1(-2 * kappa * t))
        mean = s0 * np.exp(-kappa * t)
    z = mean / np.sqrt(var)
    p = np.array([0.5 * (1 + erf(v / sqrt(2))) for v in z])
    return _binary_information(p)


def three_regimes(rng, kappa: float = 1.0, s0: float = 0.3, D: float = 0.02, n: int = 20000) -> Dict[str, object]:
    times = np.geomspace(0.05, 40.0, 18)
    rows = {}
    for label, k in (("contracting", kappa), ("flat", 0.0), ("expanding", -kappa)):
        rows[label] = {"kappa": k, "gaussian_information": gaussian_information(k, s0, D, times).tolist(),
                       "binary_exact_linear": (binary_information_linear(k, s0, D, times).tolist() if k >= 0 else None),
                       "simulation": simulate(k, s0, D, times, rng, n=n)}
    W = kappa * s0 ** 2 / D
    return {"s0": s0, "D": D, "kappa": kappa, "W": W, "plateau_gaussian": 0.5 * np.log1p(W),
            "binary_plateau": binary_plateau(W), "times": times.tolist(), "regimes": rows}


def binary_plateau(W: float) -> Dict[str, float]:
    """Accuracy and information of the sign at kappa < 0 once the expansion has outrun the noise."""
    from math import erf, sqrt
    p = 0.5 * (1 + erf(sqrt(W) / sqrt(2)))
    return {"accuracy": p, "information": float(_binary_information(np.array([p]))[0])}
