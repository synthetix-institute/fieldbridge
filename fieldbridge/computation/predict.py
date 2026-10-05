"""Predictions from the equations of a driven body, before any simulation.

The body is linearized at its steady state x* for the input offset u0: J = dF/dx, B = dF/du and H = dh/dx,
H_u = dh/du for the observables h (exact derivatives of the parsed expressions; central differences for a body given
as data). Over one interval with the input held, the deviation obeys

    s_t = Phi s_{t-1} + Gamma (u_t - u0),   Phi = e^{J dt},   Gamma = int_0^dt e^{J s} ds B,

(for a map, Phi = J and Gamma = B), and the measurement at the j-th of V times within interval t is
H (e^{J t_j} s_{t-1} + Gamma_j u_t) + H_u u_t. The impulse responses of the m V measured signals are
g_0 = [H Gamma_j + H_u]_j and g_k = [H e^{J t_j} Phi^{k-1} Gamma]_j.

- The fading class from the spectral radius of Phi over the modes that the input reaches and the observables see:
  below 1 fading memory, 1 integrating (no fading memory), above 1 unstable. Conserved amounts that the input cannot
  change are not among them.
- n_lin, the rank of [g_0, g_1, ...] at a relative tolerance of 1e-6 (the accuracy of the integration), with the
  singular values: the linear capacity at small amplitude, at most min(m V, n_O) (Gonon, Grigoryeva and Ortega 2020).
- The degree-1 profile at small amplitude, C(k) = g_k^T (sum_j g_j g_j^T)^+ g_k, summing to n_lin; for one mode
  (1 - a^2) a^(2k), a = e^{-kappa dt}.
- Odd about the steady state: F(x* - s, u0 - u) = -F(x* + s, u0 + u) and the same for the observables. With a
  symmetric input law every even degree then has zero capacity (Dambre et al. 2012).
- Linear: the drift and the observables are affine in the state and the input. Then only degree 1 has capacity.
"""
from __future__ import annotations

from typing import Dict

import numpy as np
from scipy.linalg import expm

from ..memory import spec as memory_spec
from .spec import Body

CLASSES = ("linear-memory", "odd-capacity", "nonlinear-capacity", "integrating")


def linearization(body: Body, u0: float | None = None) -> Dict:
    u0 = body.offset if u0 is None else u0
    x = body.steady(u0)
    if not hasattr(body, "real"):
        return _numeric_linearization(body, x, u0)
    sp = memory_spec._sympy()
    usym = next(s for s in body.psyms if str(s) == body.input)
    args = list(body.vsyms) + list(body.psyms)
    vals = [float(v) for v in x] + [float(body.params(u0)[str(s)]) for s in body.psyms]
    cols = list(body.vsyms) + [usym]

    def jac(exprs):
        M = sp.Matrix(exprs).jacobian(cols)
        return np.array(sp.lambdify(args, M, modules="numpy")(*vals), dtype=float).reshape(len(exprs), len(cols))

    Jd, Jh = jac(body.real.symbolic["drift"]), jac(body.obs_exprs)
    return {"x": x, "J": Jd[:, :-1], "B": Jd[:, -1], "H": Jh[:, :-1], "Hu": Jh[:, -1]}


def _numeric_linearization(body: Body, x: np.ndarray, u0: float) -> Dict:
    n = body.n
    h = 1e-6 * np.maximum(1.0, np.abs(x))
    E = np.diag(h)
    X = np.concatenate([x + E, x - E])
    U = np.full(2 * n, u0)
    fx, hx = body.f(X, U), body.observe(X, U)
    du = 1e-6 * max(1.0, abs(u0))
    xx = np.vstack([x, x])
    fu, hu = body.f(xx, np.array([u0 + du, u0 - du])), body.observe(xx, np.array([u0 + du, u0 - du]))
    return {"x": x, "J": ((fx[:n] - fx[n:]) / (2 * h[:, None])).T, "B": (fu[0] - fu[1]) / (2 * du),
            "H": ((hx[:n] - hx[n:]) / (2 * h[:, None])).T, "Hu": (hu[0] - hu[1]) / (2 * du)}


def _held(J: np.ndarray, B: np.ndarray, t: float):
    """e^{J t} and int_0^t e^{J s} ds B, from the exponential of the bordered matrix."""
    n = len(B)
    M = np.zeros((n + 1, n + 1))
    M[:n, :n], M[:n, n] = J, B
    E = expm(M * t)
    return E[:n, :n], E[:n, n]


def linear(body: Body, K: int = 400, tol: float = 1e-6) -> Dict:
    lin = linearization(body)
    J, B, H, Hu = lin["J"], lin["B"], lin["H"], lin["Hu"]
    dt, V = body.hold, body.V
    if body.form == "map":
        Phi, Gamma, nodes = J, B, [(J, B)]
    else:
        Phi, Gamma = _held(J, B, dt)
        nodes = [_held(J, B, dt * (j + 1) / V) for j in range(V)]
    g = [np.concatenate([H @ G_j + Hu for _, G_j in nodes])]
    s = Gamma.copy()
    for _ in range(K):
        g.append(np.concatenate([H @ (E_j @ s) for E_j, _ in nodes]))
        s = Phi @ s
    G = np.array(g).T
    _, sv, Vt = np.linalg.svd(G, full_matrices=False)
    n_lin = int(np.sum(sv > tol * sv[0])) if sv.size and sv[0] > 0 else 0
    profile = np.sum(Vt[:n_lin] ** 2, axis=0)
    lam, R = np.linalg.eig(Phi)
    L = np.linalg.inv(R)
    weight = np.abs(L @ Gamma) * np.linalg.norm(R, axis=0) * np.linalg.norm(H @ R, axis=0)
    live = weight > 1e-9 * weight.max() if weight.max() > 0 else np.zeros(len(lam), bool)
    rho = float(np.max(np.abs(lam[live]))) if live.any() else 0.0
    cls = "fading" if rho < 1 - 1e-9 else ("integrating" if rho < 1 + 1e-9 else "unstable")
    rates = -np.log(np.abs(lam[live]) + 1e-300) / dt if live.any() else np.array([np.inf])
    return {"class": cls, "spectral_radius": rho, "slowest_rate": float(np.min(rates)),
            "modes_reached_and_seen": int(live.sum()), "modes": len(lam), "n_lin": n_lin,
            "singular_values": (sv / sv[0]).tolist() if sv.size and sv[0] > 0 else [],
            "profile_degree_1": profile.tolist(), "total_degree_1": float(profile.sum()), "signals": body.m * V,
            "steady": lin["x"].tolist()}


def odd(body: Body, samples: int = 200, seed: int = 0) -> bool:
    """The drift and the observables are odd about the steady state, at random points around it."""
    rng = np.random.default_rng(seed)
    x0 = body.steady()
    S = rng.standard_normal((samples, body.n)) * np.maximum(np.abs(x0), 1.0) * 0.3
    U = rng.uniform(-body.amplitude, body.amplitude, samples)
    plus, minus = body.f(x0 + S, body.offset + U), body.f(x0 - S, body.offset - U)
    h0 = body.observe(x0[None, :], np.array([body.offset]))[0]
    hp = body.observe(x0 + S, body.offset + U) - h0
    hm = body.observe(x0 - S, body.offset - U) - h0
    tol = 1e-9 * (1 + np.abs(plus).max() + np.abs(hp).max())
    return bool(np.abs(plus + minus).max() < tol and np.abs(hp + hm).max() < tol)


def is_linear(body: Body, samples: int = 50, seed: int = 1) -> bool:
    """The drift and the observables are affine in the state and the input: their second differences vanish at
    random points (exact for affine maps up to rounding; a nonlinear map fails almost surely)."""
    rng = np.random.default_rng(seed)
    x0 = body.steady()
    scale = np.maximum(np.abs(x0), 1.0) * 0.3
    X = x0 + rng.standard_normal((samples, body.n)) * scale
    D = rng.standard_normal((samples, body.n)) * scale
    U = body.offset + rng.uniform(-body.amplitude, body.amplitude, samples)
    du = rng.uniform(-body.amplitude, body.amplitude, samples)
    for fn in (body.f, body.observe):
        mid = fn(X, U)
        second = fn(X + D, U + du) + fn(X - D, U - du) - 2 * mid
        if np.abs(second).max() > 1e-7 * (1 + np.abs(mid).max()):
            return False
    return True


def predict(body: Body) -> Dict:
    """The predictions and the class of the computation card."""
    out = linear(body)
    out["odd"] = odd(body)
    out["linear"] = is_linear(body)
    if out["class"] != "fading":
        out["structure"] = "integrating"
    elif out["linear"]:
        out["structure"] = "linear-memory"
    elif out["odd"]:
        out["structure"] = "odd-capacity"
    else:
        out["structure"] = "nonlinear-capacity"
    return out
