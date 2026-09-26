"""Compose memories from units on a loop, and predict before any simulation whether the loop stores,
frustrates or oscillates.

A loop of N units is a closed chain of transports: each edge maps the state of one unit to the state the edge
prefers for the next. Two transport groups cover the units here:

  sign    Z_2 = {identity, inversion}: activation or repression between genes, ferro- or antiferromagnetic
          bonds between soft spins.
  circle  affine maps of an m-fold angle: theta -> theta (alignment), theta + pi/m (anti-alignment) and
          theta -> 2 phi - theta (reflection across the bond axis, the bond-directional term of anisotropic
          capillary or dipolar rotors, phi the bond angle).

The holonomy is the composition of the transports around the loop. The loop can be satisfied by one state iff
the holonomy has a fixed point: a product of signs equal to +1; a reflection (an odd number of reflections)
always; a rotation by C only if m C = 0 mod 2 pi. The flux Phi = m C mod 2 pi measures how far it fails. This
is the frustration function of Toulouse (1977) for spins and the loop sign of Thomas (1981) for gene circuits;
for rotors the flux is set continuously by the geometry of the bonds.

Consequences that the evaluation below measures:
  reciprocal units (spins, rotors) settle into compromise states on a frustrated loop; for N equal bonds on a
  circle-group loop the lowest energy lies above the satisfied value by the fraction f = 1 - cos(Phi / N);
  non-reciprocal units (genes) cannot have two stable states on a negative loop (Thomas's rule, proved by Gouze 1998
  and Soule 2003); they relax to one state or oscillate.
"""
from __future__ import annotations

import time
from typing import Dict, List, Optional, Sequence

import numpy as np

from . import analysis as an
from .carriers import euclid, orthant
from .identity import Realization
from .library import rotor_graph


# ------------------------------------------------------------------------------------------------ holonomy
def holonomy(loop: Dict) -> Dict[str, object]:
    if loop["group"] == "sign":
        S = int(np.prod(loop["signs"]))
        return {"group": "sign", "holonomy": "+1" if S > 0 else "-1", "satisfiable": S > 0,
                "flux": 0.0 if S > 0 else float(np.pi)}
    m = loop["m"]
    S, C = 1, 0.0
    for kind, phi in zip(loop["kinds"], loop["phis"]):
        s, c = {"align": (1, 0.0), "anti": (1, np.pi / m), "reflect": (-1, 2.0 * phi)}[kind]
        S, C = s * S, s * C + c
    if S < 0:
        return {"group": "circle", "holonomy": "reflection", "satisfiable": True, "flux": 0.0}
    flux = float(np.angle(np.exp(1j * m * C)))
    return {"group": "circle", "holonomy": "rotation", "satisfiable": bool(abs(flux) < 1e-9), "flux": flux}


def predicted_frustration(loop: Dict) -> float:
    """Excess energy fraction of the best compromise on a circle-group loop of equal bonds."""
    hol = holonomy(loop)
    n = len(loop["kinds"])
    return 0.0 if hol["satisfiable"] else float(1.0 - np.cos(hol["flux"] / n))


# ------------------------------------------------------------------------------------------------ loops
def _ring_positions(N: int, rng=None, jitter: float = 0.0) -> np.ndarray:
    a = 2 * np.pi * np.arange(N) / N
    pos = np.stack([np.cos(a), np.sin(a)], axis=1) / (2 * np.sin(np.pi / N))
    if rng is not None and jitter > 0:
        pos = pos + rng.normal(scale=jitter, size=pos.shape)
    return pos


def gene_loop(N: int, repress: Sequence[bool], alpha: float = 10.0, n: float = 4.0, noise: float = 0.05) -> Realization:
    """du_i/dt = alpha g_i(u_{i-1}) - u_i, with g = 1/(1 + u^n) for repression and u^n/(1 + u^n) for activation."""
    rep = np.asarray(repress, bool)

    def drift(q, p):
        x = np.roll(q, 1, axis=-1) ** p["n"]
        return p["alpha"] * np.where(rep, 1.0 / (1.0 + x), x / (1.0 + x)) - q

    real = Realization(
        name=f"gene loop, N = {N}, {int(rep.sum())} repression(s)",
        carrier=orthant(f"protein concentrations u_1..u_{N}", N, scale=alpha), drift=drift,
        params={"alpha": alpha, "n": n}, control="alpha", control_range=(0.5, 12.0), noise=noise, dt=0.01,
        closure="dilution; each gene regulated by its predecessor", observable="which genes are expressed",
        provenance="transcriptional feedback loop", transport="sign: activation +1, repression -1",
        tags=("biology", "non-reciprocal", "loop"))
    real.loop = {"group": "sign", "signs": [-1 if r else 1 for r in rep], "reciprocal": False}
    real.geometry = {"pos": _ring_positions(N).tolist(), "src": list(range(N)), "tgt": [(i + 1) % N for i in range(N)],
                     "edge_kind": ["repress" if r else "activate" for r in np.roll(rep, -1)]}
    return real


def spin_loop(N: int, signs: Sequence[int], k: float = 0.4, eps: float = 1.0, noise: float = 0.05) -> Realization:
    """Soft spins, dq_i/dt = eps q_i - q_i^3 + k (s_i q_{i+1} + s_{i-1} q_{i-1}), s_i = +1 (ferro) or -1 (antiferro)."""
    sg = np.asarray(signs, float)

    def drift(q, p):
        return (p["eps"] * q - q ** 3
                + p["k"] * (sg * np.roll(q, -1, axis=-1) + np.roll(sg, 1) * np.roll(q, 1, axis=-1)))

    def potential(q, p):
        return np.sum(-0.5 * p["eps"] * q ** 2 + 0.25 * q ** 4 - p["k"] * sg * q * np.roll(q, -1, axis=-1), axis=-1)

    real = Realization(
        name=f"spin loop, N = {N}, {int((sg < 0).sum())} antiferromagnetic bond(s)",
        carrier=euclid(f"soft spins q_1..q_{N}", N, scale=1.5), drift=drift, potential=potential,
        params={"k": k, "eps": eps}, control="eps", control_range=(-0.5, 1.5), noise=noise, dt=0.01,
        closure="thermal bath", observable="pattern of spin signs", provenance="Ising-like soft spins",
        transport="sign: ferro +1, antiferro -1", tags=("gradient", "reciprocal", "loop"))
    real.loop = {"group": "sign", "signs": [int(x) for x in sg], "reciprocal": True}
    real.geometry = {"pos": _ring_positions(N).tolist(), "src": list(range(N)), "tgt": [(i + 1) % N for i in range(N)],
                     "edge_kind": ["ferro" if x > 0 else "antiferro" for x in sg]}
    return real


def rotor_loop(N: int, kinds: Sequence[str], positions: np.ndarray, m: int = 2, J: float = 1.0,
               noise: float = 0.2) -> Realization:
    """Rotors at fixed positions on a loop; each bond is alignment, anti-alignment or bond-directional."""
    pos = np.asarray(positions, float)
    src = np.arange(N)
    tgt = (src + 1) % N
    d = pos[tgt] - pos[src]
    phi = np.arctan2(d[:, 1], d[:, 0])
    wa = np.array([{"align": 1.0, "anti": -1.0, "reflect": 0.0}[k] for k in kinds])
    wc = np.array([1.0 if k == "reflect" else 0.0 for k in kinds])
    real = rotor_graph(f"rotor loop, N = {N}, m = {m}: " + " ".join(k[:3] for k in kinds), m, src, tgt, phi, wa, wc,
                       J, J, 1.0, noise, dt=0.01,
                       provenance="caged anisotropic rotors (capillary for m = 2, dipolar for m = 1)")
    real.control = None
    real.loop = {"group": "circle", "m": m, "kinds": list(kinds), "phis": phi.tolist(), "reciprocal": True}
    real.geometry = {"pos": pos.tolist(), "src": src.tolist(), "tgt": tgt.tolist(), "edge_kind": list(kinds)}
    # alignment keeps theta_i - theta_j and reflection keeps theta_i + theta_j: with an even number of
    # reflections the rotation theta_i -> theta_i + e_i a (e_i = +-1, flipping at each reflection) is exact
    signs = np.cumprod([1.0] + [(-1.0 if k == "reflect" else 1.0) for k in kinds[:-1]])
    if int(np.sum([k == "reflect" for k in kinds])) % 2 == 0:
        period = 2 * np.pi / m
        real.null_modes = signs[None, :]
        real.canonical = lambda q, e=signs: np.mod(q - e * (q[..., :1] * e[0]), period)
    return real


# ------------------------------------------------------------------------------------------------ evaluation
def evaluate_loop(real: Realization, rng, n_starts: int = 32) -> Dict[str, object]:
    t0 = time.time()
    loop = real.loop
    hol = holonomy(loop)
    states, counts, spectra, n_unconv, continuum = an.stored_states(real, rng, n_starts)
    row: Dict[str, object] = {"name": real.name, "group": loop["group"], "reciprocal": loop["reciprocal"],
                              "N": real.carrier.dim, "holonomy": hol, "stable_states": len(states),
                              "unsettled_starts": n_unconv, "geometry": real.geometry}
    if loop["group"] == "sign":
        row["signs"] = list(loop["signs"])
    if getattr(real, "showcase", None):
        row["showcase"] = real.showcase
    if loop["group"] == "sign" and not loop["reciprocal"]:
        row["prediction"] = "at most one stable state" if not hol["satisfiable"] else "two or more stable states possible"
        row["measured"] = ("oscillation" if len(states) == 0 and n_unconv > 0 else f"{len(states)} stable state(s)")
        row["consistent"] = bool(hol["satisfiable"] or len(states) <= 1)
        row["stores"] = len(states) >= 2
    elif loop["group"] == "sign":
        E = real.V(np.array(states))
        g = states[int(np.argmin(E))]
        sg = np.asarray(loop["signs"], float)
        unsat = int(np.sum(sg * g * np.roll(g, -1) < 0))
        row.update(prediction="frustrated" if not hol["satisfiable"] else "satisfied", unsatisfied_bonds=unsat,
                   ground_state=g.tolist(), measured="frustrated" if unsat else "satisfied",
                   consistent=bool((unsat > 0) == (not hol["satisfiable"])), stores=len(states) >= 2,
                   ground_degeneracy=int(np.sum(E < E.min() + 1e-6)))
    else:
        E = real.V(np.array(states))
        g = states[int(np.argmin(E))]
        wsum = float(np.sum(np.abs(real.graph["wa"]) + np.abs(real.graph["wc"])))
        f = float((E.min() + wsum) / wsum)
        fp = predicted_frustration(loop)
        row.update(prediction=fp, measured=f, frustration=f, consistent=bool(abs(f - fp) < 0.01),
                   ground_state=g.tolist(), stores=len(states) >= 2,
                   continuous_symmetry=real.null_modes is not None)
    row["seconds"] = time.time() - t0
    return row


def showcase_loops(rng) -> List[Realization]:
    """Three rotor loops for drawings: a triangle (odd: a reflection, satisfied), a square (even, mismatch 0) and a
    quadrilateral whose loop mismatch is close to pi (the most frustrated four-loop)."""
    tri = rotor_loop(3, ["reflect"] * 3, _ring_positions(3))
    sq = rotor_loop(4, ["reflect"] * 4, _ring_positions(4))
    best, quad = None, None
    for _ in range(400):
        pos = _ring_positions(4, rng, jitter=0.35)
        cand = rotor_loop(4, ["reflect"] * 4, pos)
        gap = abs(abs(holonomy(cand.loop)["flux"]) - np.pi)
        if best is None or gap < best:
            best, quad = gap, cand
        if gap < 0.02:
            break
    for r, label in ((tri, "triangle: odd loop, satisfied"), (sq, "square: mismatch 0"),
                     (quad, "quadrilateral: mismatch near pi")):
        r.showcase = label
    return [tri, sq, quad]


def standard_loops(rng, rotor_shapes: int = 8, mixed: int = 16) -> List[Realization]:
    loops: List[Realization] = showcase_loops(rng)
    for N in (2, 3, 4, 5):
        for r in range(N + 1):
            loops.append(gene_loop(N, [i < r for i in range(N)]))
    for N in (3, 4, 5, 6):
        for r in range(N + 1):
            loops.append(spin_loop(N, [-1 if i < r else 1 for i in range(N)]))
    for N in (3, 4, 5, 6):
        for _ in range(rotor_shapes):
            loops.append(rotor_loop(N, ["reflect"] * N, _ring_positions(N, rng, jitter=0.3)))
    for _ in range(mixed):
        N = int(rng.integers(3, 7))
        kinds = list(rng.choice(["align", "anti", "reflect"], size=N))
        loops.append(rotor_loop(N, kinds, _ring_positions(N, rng, jitter=0.3)))
    return loops


def run_loops(rng, loops: Optional[List[Realization]] = None) -> List[Dict[str, object]]:
    loops = loops if loops is not None else standard_loops(rng)
    return [evaluate_loop(real, rng) for real in loops]
