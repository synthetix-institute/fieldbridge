"""A library of realizations on different carriers, each with its physical provenance.

Every entry is a Realization: a drift on a carrier with material parameters, a control parameter that a
protocol can sweep, a default noise and, for gradient systems, the potential. The entries span:

  pitchfork          R^1      Landau normal form; the reference write (at the ridge) and retention (in the wells)
  bistable_ring      R^n      bistable units coupled by a discrete Laplacian (sequence-memory substrate)
  rotor_patch        (S^1)^N  caged anisotropic colloids at a fluid interface (capillary quadrupoles), or
                              in-plane dipoles at the same positions
  dipole_exact_map   (S^1)^N  the dipole model with doubled bond angles, the exact image of rotor_patch
  schlogl            R_+      autocatalytic chemistry with two stable concentrations (mass action)
  adaptive_network   R_+^K    conductances of K parallel tubes that adapt to the flux they carry (Physarum)
  toggle             R_+^2    two mutually repressing genes (Gardner, Cantor & Collins 2000)
  repressor_ring     R_+^N    N genes in a ring of repression (N = 3: repressilator)
  compartments       R^2      Hyperion's three-compartment exchange: an observable x and an unobserved imbalance z
"""
from __future__ import annotations

import math
from typing import Dict, Tuple

import numpy as np

from .carriers import euclid, orthant, torus
from .identity import Realization


# --------------------------------------------------------------------------------------------------
def pitchfork(eps: float = 1.0, noise: float = 0.08) -> Realization:
    return Realization(
        name="pitchfork, dx/dt = eps x - x^3",
        carrier=euclid("order parameter x", 1, scale=2.0),
        drift=lambda q, p: p["eps"] * q - q ** 3,
        potential=lambda q, p: np.sum(-0.5 * p["eps"] * q ** 2 + 0.25 * q ** 4, axis=-1),
        params={"eps": eps}, control="eps", control_range=(-0.5, 1.5), noise=noise, dt=0.01,
        closure="overdamped coordinate in a thermal bath",
        observable="sign of x", provenance="Landau normal form of a symmetry-breaking transition",
        tags=("gradient", "normal form"))


def bistable_ring(n: int = 8, k: float = 0.3, eps: float = 1.0, noise: float = 0.06) -> Realization:
    def drift(q, p):
        lap = np.roll(q, 1, axis=-1) + np.roll(q, -1, axis=-1) - 2.0 * q
        return p["eps"] * q - q ** 3 + p["k"] * lap

    def potential(q, p):
        return (np.sum(-0.5 * p["eps"] * q ** 2 + 0.25 * q ** 4, axis=-1)
                + 0.5 * p["k"] * np.sum((np.roll(q, -1, axis=-1) - q) ** 2, axis=-1))

    return Realization(
        name=f"ring of {n} bistable coordinates, Laplacian coupling k",
        carrier=euclid(f"coordinates q_1..q_{n}", n, scale=1.5), drift=drift, potential=potential,
        params={"eps": eps, "k": k}, control="eps", control_range=(-0.5, 1.5), noise=noise, dt=0.01,
        closure="ring (periodic), thermal bath", observable="pattern of signs",
        provenance="minimal substrate of order memory (three-write criterion)", transport="identity (ferromagnetic)",
        tags=("gradient",))


# --------------------------------------------------------------------------------------------------
def _patch(rows: int = 3, cols: int = 4, disorder: float = 0.16, cutoff: float = 1.45, seed: int = 17):
    rng = np.random.default_rng(seed)
    pos = np.array([[c + 0.5 * (r % 2), r * math.sqrt(3) / 2] for r in range(rows) for c in range(cols)], float)
    pos += rng.normal(scale=disorder, size=pos.shape)
    src, tgt = np.triu_indices(len(pos), 1)
    d = pos[tgt] - pos[src]
    dist = np.linalg.norm(d, axis=1)
    keep = dist < cutoff
    src, tgt, d, dist = src[keep], tgt[keep], d[keep], dist[keep]
    return pos, src, tgt, dist, np.arctan2(d[:, 1], d[:, 0])


def rotor_graph(name: str, m: int, src, tgt, phi, wa, wc, J0: float, G0: float, lam: float, noise: float,
                mobility: float = 1.0, dt: float = 0.002, provenance: str = "", closure: str = "",
                extra: Dict[str, float] | None = None) -> Realization:
    """U = -sum_edges [J wa cos(m (q_i - q_j)) + G wc cos(m (q_i + q_j - 2 phi))], J = J0 lam, G = G0 lam."""
    n = int(max(src.max(), tgt.max())) + 1

    def forces(q, p):
        J, G, mu = p["J0"] * p["lam"], p["G0"] * p["lam"], p["mobility"]
        d = q[:, src] - q[:, tgt]
        s = q[:, src] + q[:, tgt] - 2.0 * phi
        a = J * wa * m * np.sin(m * d)
        c = G * wc * m * np.sin(m * s)
        out = np.zeros_like(q)
        np.add.at(out, (slice(None), src), -(a + c))
        np.add.at(out, (slice(None), tgt), a - c)
        return mu * out

    def potential(q, p):
        J, G = p["J0"] * p["lam"], p["G0"] * p["lam"]
        d = q[:, src] - q[:, tgt]
        s = q[:, src] + q[:, tgt] - 2.0 * phi
        return -np.sum(J * wa * np.cos(m * d) + G * wc * np.cos(m * s), axis=-1)

    params = {"J0": J0, "G0": G0, "lam": lam, "mobility": mobility, **(extra or {})}
    real = Realization(
        name=name, carrier=torus(f"orientations of {n} rotors", n, 2.0 * np.pi / m), drift=forces,
        potential=potential, params=params, control="lam", control_range=(0.05, 1.5), noise=noise, dt=dt,
        closure=closure or "positions caged; rotational Brownian motion", observable="overlap with a stored angular pattern",
        provenance=provenance, transport="reflection across the bond axis (bond-directional term) + identity (alignment)",
        tags=("gradient", "network"))
    real.graph = {"src": src, "tgt": tgt, "phi": phi, "wa": wa, "wc": wc, "m": m}
    return real


def rotor_patch(kind: str = "capillary", lam: float = 0.6, noise: float = 1.0, seed: int = 17) -> Realization:
    """Twelve caged rotors on a jittered triangular patch.

    capillary: nematic rods at a fluid interface, U = -sum [J cos 2(dtheta) + g (r0/r)^4 cos 2(sum - 2 phi)],
               J : g = 4 : 5 (the colloid Letter); period pi.
    dipolar:   in-plane point dipoles at the same positions, E = (mu^2/r^3)[cos(dpsi) - 3 cos(psi_i - phi)cos(psi_j - phi)]
               = -(mu^2/r^3)[1/2 cos(dpsi) + 3/2 cos(psi_i + psi_j - 2 phi)]; period 2 pi.
    """
    pos, src, tgt, dist, phi = _patch(seed=seed)
    r0 = float(np.median(dist))
    if kind == "capillary":
        real = rotor_graph("caged capillary rotors (nematic, m = 2)", 2, src, tgt, phi, np.ones_like(dist),
                           np.clip((r0 / dist) ** 4, 0, 6), 4.0, 5.0, lam, noise,
                           provenance="anisotropic colloids at a fluid interface, quadrupolar capillary coupling")
    elif kind == "dipolar":
        w = (r0 / dist) ** 3
        real = rotor_graph("caged in-plane dipoles (vector, m = 1)", 1, src, tgt, phi, w, w, 0.5, 1.5, lam, noise,
                           provenance="magnetic or electric point dipoles with fixed centres")
    else:
        raise ValueError(kind)
    real.geometry = {"pos": pos.tolist(), "src": src.tolist(), "tgt": tgt.tolist()}
    return real


def dipole_exact_map(source: Realization) -> Realization:
    """Image of a nematic rotor graph (m = 2) under psi = 2 theta: a vector model with bond angles 2 phi."""
    g = source.graph
    real = rotor_graph("vector rotors with doubled bond angles (exact image of the nematic network)", 1, g["src"], g["tgt"],
                       2.0 * g["phi"], g["wa"], g["wc"], source.params["J0"], source.params["G0"], source.params["lam"],
                       4.0 * source.noise, mobility=4.0 * source.params["mobility"], dt=source.dt,
                       provenance="abstract carrier: exact map of the capillary network (not realizable by positions)")
    real.geometry = source.geometry
    return real


# --------------------------------------------------------------------------------------------------
def schlogl(b: float = 1.875, a: float = 4.5, k3: float = 5.75, noise: float = 0.04) -> Realization:
    """Schlogl's second model, A + 2X <-> 3X, X <-> B, with unit rate constants for the cubic step:
    dx/dt = -x^3 + a x^2 - k3 x + b, where a and b are set by the bath concentrations of A and B.
    At the default values the stable concentrations are 0.5 and 2.5."""
    return Realization(
        name="autocatalytic chemistry (Schlogl), dx/dt = -x^3 + a x^2 - k3 x + b",
        carrier=orthant("concentration x", 1, scale=3.0),
        drift=lambda q, p: -q ** 3 + p["a"] * q ** 2 - p["k3"] * q + p["b"],
        potential=lambda q, p: np.sum(q ** 4 / 4 - p["a"] * q ** 3 / 3 + p["k3"] * q ** 2 / 2 - p["b"] * q, axis=-1),
        params={"b": b, "a": a, "k3": k3}, control="b", control_range=(1.075, 2.675), noise=noise, dt=0.005,
        closure="well-mixed reactor fed by baths of A and B", observable="concentration x",
        provenance="A + 2X <-> 3X, X <-> B (Schlogl 1972)", material=("a",), tags=("gradient", "chemistry"))


def adaptive_network(K: int = 2, mu: float = 2.0, lengths=None, noise: float = 0.005) -> Realization:
    """K parallel tubes carry a fixed total flux between a source and a sink; each conductance D_k grows with
    the flux it carries and decays otherwise, dD_k/dt = f(|Q_k|) - D_k with f(Q) = Q^mu / (1 + Q^mu) and
    Q_k = (D_k/L_k) / sum_j (D_j/L_j) (Kirchhoff), after Tero, Kobayashi & Nakagaki (2007)."""
    L = np.ones(K) if lengths is None else np.asarray(lengths, float)

    def drift(q, p):
        Ls = np.array([p[f"L{k + 1}"] for k in range(K)])
        g = np.abs(q) / Ls
        Q = g / np.maximum(g.sum(axis=-1, keepdims=True), 1e-12)
        fQ = Q ** p["mu"]
        return fQ / (1.0 + fQ) - q

    return Realization(
        name=f"adaptive transport network, {K} parallel tubes",
        carrier=orthant(f"tube conductances D_1..D_{K}", K, scale=1.0), drift=drift,
        params={"mu": mu, **{f"L{k + 1}": float(L[k]) for k in range(K)}}, control="mu", control_range=(0.5, 3.0),
        noise=noise, dt=0.02, closure="fixed total flux between a source and a sink (Kirchhoff's laws)",
        observable="which tube carries the flow",
        provenance="Physarum-type adaptation of tube conductance to flux (Tero, Kobayashi & Nakagaki 2007)",
        material=tuple(f"L{k + 1}" for k in range(1, K)), tags=("network", "parameter carrier"))


def toggle(alpha: float = 4.0, n: float = 2.0, noise: float = 0.05) -> Realization:
    def drift(q, p):
        u, v = q[:, 0], q[:, 1]
        return np.stack([p["alpha"] / (1 + v ** p["n"]) - u, p["alpha"] / (1 + u ** p["n"]) - v], axis=-1)

    return Realization(
        name="genetic toggle, du/dt = a/(1+v^n) - u, dv/dt = a/(1+u^n) - v",
        carrier=orthant("repressor concentrations (u, v)", 2, scale=alpha), drift=drift,
        params={"alpha": alpha, "n": n}, control="alpha", control_range=(0.5, 6.0), noise=noise, dt=0.01,
        closure="cell with dilution; intrinsic noise", observable="which repressor is high",
        provenance="mutual repression of two genes (Gardner, Cantor & Collins 2000)", transport="sign inversion (repression)",
        tags=("biology", "non-gradient"))


def repressor_ring(N: int = 3, alpha: float = 10.0, n: float = 4.0, noise: float = 0.05) -> Realization:
    def drift(q, p):
        return p["alpha"] / (1 + np.roll(q, 1, axis=-1) ** p["n"]) - q

    return Realization(
        name=f"ring of {N} repressors, du_i/dt = a/(1+u_(i-1)^n) - u_i",
        carrier=orthant(f"repressor concentrations u_1..u_{N}", N, scale=alpha), drift=drift,
        params={"alpha": alpha, "n": n}, control="alpha", control_range=(0.5, 12.0), noise=noise, dt=0.01,
        closure="cell with dilution; intrinsic noise", observable="pattern of high repressors",
        provenance="repressilator (Elowitz & Leibler 2000) for N = 3", transport="sign inversion (repression)",
        tags=("biology", "non-gradient"))


def compartments(ax: float = 7.5, az: float = 3.9, b: float = 1.5, c: float = 4.5, noise: float = 0.02) -> Realization:
    """Hyperion's exchange model: x = p0 - 1/3 observed, z = p1 - p2 unobserved; dx = -ax x + b z, dz = c x - az z."""
    def drift(q, p):
        x, z = q[:, 0], q[:, 1]
        return np.stack([-p["ax"] * x + p["b"] * z, p["c"] * x - p["az"] * z], axis=-1)

    real = Realization(
        name="three compartments: observable x, unobserved imbalance z",
        carrier=euclid("compartment imbalances (x, z)", 2, scale=0.3), drift=drift,
        params={"ax": ax, "az": az, "b": b, "c": c}, control=None, noise=noise, dt=0.002,
        closure="probability conserved; symmetric exchange", observable="return x to compartment 0 (z not observed)",
        observable_fn=lambda q: q[..., 0], provenance="Hyperion, three-compartment memory kernel",
        tags=("linear", "unobserved block"))
    real.linear = {"A": np.array([[-ax]]), "B": np.array([[b]]), "C": np.array([[c]]), "D": np.array([[-az]])}
    return real


def catalog() -> Tuple[Realization, ...]:
    cap = rotor_patch("capillary")
    return (pitchfork(), bistable_ring(), cap, rotor_patch("dipolar"), dipole_exact_map(cap), schlogl(),
            adaptive_network(), adaptive_network(lengths=(1.0, 1.15)), toggle(), repressor_ring(3), repressor_ring(4),
            compartments())
