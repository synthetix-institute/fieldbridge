"""Networks of any topology: predictions from structure alone, and their check by simulation.

A network is a set of units joined by edges, each edge carrying a transport, the map that takes the state of one
unit to the state the edge prefers for the other:

  rotor      m-fold angles; alignment theta_j = theta_i, anti-alignment theta_j = theta_i + pi/m, reflection across
             the bond axis theta_j = 2 phi_ij - theta_i (bond-directional coupling). Reciprocal (gradient).
  soft spin  bistable coordinates; ferromagnetic (+1) or antiferromagnetic (-1) bonds. Reciprocal.
  gene       concentrations; each gene is produced at a rate set by its regulators, activation (+1) or
             repression (-1). Directed, not reciprocal.

Predictions that use only the edges and their transports (no integration of the dynamics):

  rotors     the network can satisfy every bond iff the transports admit a consistent assignment. Propagating
             theta_i = S_i theta_0 + C_i along a spanning tree turns every other edge into a constraint on theta_0:
             a rotation-type edge (S_j = s S_i) requires a vanishing flux, a reflection-type edge fixes theta_0 up
             to half a period, and all such values must agree. If no constraint is of reflection type, the signed
             rotation theta_i -> theta_i + S_i a leaves every bond unchanged: a continuous symmetry, a flat
             direction and Law-2 holding.
  spins      satisfiable iff the signed graph is balanced (every cycle positive; Harary 1953), which a two-colouring
             search decides.
  genes      without a positive directed cycle there is at most one steady state (Thomas's rule; Soule 2003); if
             every undirected cycle of the interaction graph is positive the system is monotone and its
             trajectories converge to steady states (Hirsch; Smith 1995), so it cannot oscillate.
"""
from __future__ import annotations

import itertools
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import analysis as an
from .carriers import euclid, orthant
from .identity import Realization
from .library import rotor_graph

TRANSPORT = {"align": (1, 0.0), "anti": (1, None), "reflect": (-1, None)}


# ------------------------------------------------------------------------------------------------ structure
def rotor_constraints(n: int, edges: Sequence[Tuple[int, int, str]], phis: Sequence[float], m: int) -> Dict[str, object]:
    """Exact satisfiability of a rotor network from its transports (no dynamics)."""
    P = 2.0 * np.pi / m
    adj: Dict[int, List[Tuple[int, int]]] = {i: [] for i in range(n)}
    for k, (i, j, _) in enumerate(edges):
        adj[i].append((j, k))
        adj[j].append((i, k))

    def edge_map(k: int, frm: int) -> Tuple[int, float]:
        i, j, kind = edges[k]
        if kind == "align":
            return 1, 0.0
        if kind == "anti":
            return 1, np.pi / m
        return -1, 2.0 * phis[k]  # reflection is its own inverse: the same map in both directions

    S = np.zeros(n, int)
    C = np.zeros(n)
    root_of = -np.ones(n, int)
    tree = set()
    for r in range(n):
        if root_of[r] >= 0:
            continue
        root_of[r], S[r], C[r] = r, 1, 0.0
        stack = [r]
        while stack:
            i = stack.pop()
            for j, k in adj[i]:
                if root_of[j] < 0:
                    s, c = edge_map(k, i)
                    root_of[j], S[j], C[j] = r, s * S[i], s * C[i] + c
                    tree.add(k)
                    stack.append(j)
    fluxes, pins, frustrated = [], {}, False
    for k, (i, j, kind) in enumerate(edges):
        if k in tree:
            continue
        s, c = edge_map(k, i)
        coef = S[j] - s * S[i]
        rhs = s * C[i] + c - C[j]
        if coef == 0:
            flux = float(np.angle(np.exp(1j * m * rhs)))
            fluxes.append({"edge": k, "flux": flux})
            if abs(flux) > 1e-9:
                frustrated = True
        else:
            x = (rhs / coef) % (P / 2.0)
            pins.setdefault(int(root_of[i]), []).append(float(x))
    for r, xs in pins.items():
        xs = np.asarray(xs)
        d = np.abs((xs[:, None] - xs[None, :] + P / 4.0) % (P / 2.0) - P / 4.0)
        if np.any(d > 1e-9):
            frustrated = True
    roots = sorted(set(int(x) for x in root_of))
    free = [r for r in roots if r not in pins]
    return {"satisfiable": not frustrated, "rotation_fluxes": fluxes,
            "continuous_symmetry": bool(free), "signed_rotation": [int(s) for s in S],
            "components": len(roots)}


def sign_balance(n: int, edges: Sequence[Tuple[int, int, int]]) -> Dict[str, object]:
    """Harary balance of an undirected signed graph by two-colouring."""
    adj: Dict[int, List[Tuple[int, int]]] = {i: [] for i in range(n)}
    for i, j, s in edges:
        adj[i].append((j, s))
        adj[j].append((i, s))
    colour = [0] * n
    balanced = True
    for r in range(n):
        if colour[r]:
            continue
        colour[r] = 1
        stack = [r]
        while stack:
            i = stack.pop()
            for j, s in adj[i]:
                want = colour[i] * s
                if colour[j] == 0:
                    colour[j] = want
                    stack.append(j)
                elif colour[j] != want:
                    balanced = False
    return {"balanced": balanced, "gauge": colour}


def directed_cycles(n: int, edges: Sequence[Tuple[int, int, int]], max_len: int = 10) -> List[Tuple[List[int], int]]:
    """Simple directed cycles with their signs (product of edge signs), for small graphs."""
    out_edges: Dict[int, List[Tuple[int, int]]] = {i: [] for i in range(n)}
    for i, j, s in edges:
        out_edges[i].append((j, s))
    cycles = []
    for start in range(n):
        stack = [(start, [start], 1)]
        while stack:
            node, path, sign = stack.pop()
            for nxt, s in out_edges[node]:
                if nxt == start:
                    cycles.append((list(path), sign * s))
                elif nxt > start and nxt not in path and len(path) < max_len:
                    stack.append((nxt, path + [nxt], sign * s))
    return cycles


# ------------------------------------------------------------------------------------------------ builders
def _edge_terms(edge) -> Dict[str, float]:
    """An edge is (i, j, kind) with unit strength, or (i, j, {"align": a, "anti": b, "reflect": g})."""
    i, j, spec = edge
    if isinstance(spec, str):
        return {spec: 1.0}
    return {k: float(v) for k, v in spec.items() if float(v) != 0.0}


def pure_kind(edge) -> Optional[str]:
    terms = _edge_terms(edge)
    return next(iter(terms)) if len(terms) == 1 else None


def rotor_network(positions: np.ndarray, edges: Sequence, m: int = 2, J: float = 1.0,
                  noise: float = 0.2) -> Realization:
    pos = np.asarray(positions, float)
    src = np.array([e[0] for e in edges])
    tgt = np.array([e[1] for e in edges])
    d = pos[tgt] - pos[src]
    phi = np.arctan2(d[:, 1], d[:, 0])
    terms = [_edge_terms(e) for e in edges]
    wa = np.array([t.get("align", 0.0) - t.get("anti", 0.0) for t in terms])
    wc = np.array([t.get("reflect", 0.0) for t in terms])
    real = rotor_graph(f"rotor network, {len(pos)} rotors, {len(edges)} bonds", m, src, tgt, phi, wa, wc, J, J, 1.0,
                       noise, dt=0.01, provenance="caged anisotropic rotors with mixed bonds")
    real.control = None
    kinds = [pure_kind(e) or "mixed" for e in edges]
    real.geometry = {"pos": pos.tolist(), "src": src.tolist(), "tgt": tgt.tolist(), "edge_kind": kinds}
    real.network = {"type": "rotor", "m": m, "edges": [[int(e[0]), int(e[1]), e[2]] for e in edges],
                    "phis": phi.tolist(), "pure": all(k != "mixed" for k in kinds)}
    if real.network["pure"]:
        cons = rotor_constraints(len(pos), [(int(e[0]), int(e[1]), pure_kind(e)) for e in edges], phi, m)
        if cons["continuous_symmetry"]:
            signs = np.asarray(cons["signed_rotation"], float)
            period = 2 * np.pi / m
            real.null_modes = signs[None, :]
            real.canonical = lambda q, e=signs: np.mod(q - e * (q[..., :1] * e[0]), period)
    return real


def spin_network(n: int, edges: Sequence[Tuple[int, int, int]], k: float = 0.4, eps: float = 1.0,
                 noise: float = 0.05) -> Realization:
    src = np.array([e[0] for e in edges])
    tgt = np.array([e[1] for e in edges])
    sg = np.array([e[2] for e in edges], float)

    def drift(q, p):
        out = p["eps"] * q - q ** 3
        pull = np.zeros_like(q)
        np.add.at(pull, (slice(None), src), p["k"] * sg * q[:, tgt])
        np.add.at(pull, (slice(None), tgt), p["k"] * sg * q[:, src])
        return out + pull

    def potential(q, p):
        return (np.sum(-0.5 * p["eps"] * q ** 2 + 0.25 * q ** 4, axis=-1)
                - p["k"] * np.sum(sg * q[:, src] * q[:, tgt], axis=-1))

    real = Realization(name=f"soft-spin network, {n} spins, {len(edges)} bonds",
                       carrier=euclid(f"soft spins q_1..q_{n}", n, scale=1.5), drift=drift, potential=potential,
                       params={"k": k, "eps": eps}, control="eps", control_range=(-0.5, 1.5), noise=noise, dt=0.01,
                       closure="thermal bath", observable="pattern of spin signs",
                       transport="sign: ferro +1, antiferro -1", tags=("gradient", "network"))
    real.network = {"type": "spin", "edges": [list(e) for e in edges]}
    return real


def gene_network(n: int, regulations: Sequence[Tuple[int, int, int]], alpha: float = 10.0, hill: float = 4.0,
                 noise: float = 0.05) -> Realization:
    """du_j/dt = alpha prod_{i -> j} h_s(u_i) - u_j, with h_+ = u^n/(1+u^n) and h_- = 1/(1+u^n)."""
    regs = [(int(i), int(j), int(s)) for i, j, s in regulations]

    def drift(q, p):
        prod = np.ones_like(q)
        for i, j, s in regs:
            x = q[:, i] ** p["n"]
            prod[:, j] *= (x / (1.0 + x)) if s > 0 else (1.0 / (1.0 + x))
        return p["alpha"] * prod - q

    real = Realization(name=f"gene network, {n} genes, {len(regs)} regulations",
                       carrier=orthant(f"protein concentrations u_1..u_{n}", n, scale=alpha), drift=drift,
                       params={"alpha": alpha, "n": hill}, control="alpha", control_range=(0.5, 12.0), noise=noise,
                       dt=0.01, closure="dilution; production set by the regulators", observable="which genes are on",
                       transport="sign: activation +1, repression -1", tags=("biology", "non-reciprocal", "network"))
    real.network = {"type": "gene", "edges": [list(e) for e in regs]}
    return real


# ------------------------------------------------------------------------------------------------ predictions
def predict(real: Realization) -> Dict[str, object]:
    """What the structure alone predicts; every entry names the structural input it used."""
    net = real.network
    n = real.carrier.dim
    if net["type"] == "rotor":
        if not net.get("pure", True):
            return {"uses": "bond transports", "satisfiable": None,
                    "note": "bonds combine several transports; satisfiability is not decided by the transports alone"}
        c = rotor_constraints(n, [(int(e[0]), int(e[1]), pure_kind(e)) for e in net["edges"]], net["phis"], net["m"])
        return {"uses": "bond transports and bond angles", "satisfiable": c["satisfiable"],
                "continuous_symmetry": c["continuous_symmetry"],
                "loss_law": ("diffusion along the signed rotation (Law 2)" if c["continuous_symmetry"]
                             else "activation between discrete states (Law 3)"),
                "rotation_fluxes": c["rotation_fluxes"]}
    if net["type"] == "spin":
        b = sign_balance(n, [tuple(e) for e in net["edges"]])
        return {"uses": "bond signs", "satisfiable": b["balanced"],
                "loss_law": "activation between discrete states (Law 3)"}
    edges = [tuple(e) for e in net["edges"]]
    cyc = directed_cycles(n, edges)
    positive = any(s > 0 for _, s in cyc)
    negative = any(s < 0 for _, s in cyc)
    monotone = sign_balance(n, edges)["balanced"]
    if not positive:
        states = "at most one steady state (no positive loop)"
    else:
        states = "two or more stable states possible (a positive loop)"
    return {"uses": "signs of the regulations", "positive_loop": positive, "negative_loop": negative,
            "monotone": monotone, "states": states,
            "oscillation": ("impossible (monotone: every undirected cycle positive)" if monotone
                            else ("possible (a negative loop)" if negative else "not excluded"))}


# ------------------------------------------------------------------------------------------------ measurement
def measure(real: Realization, rng, n_starts: int = 32) -> Dict[str, object]:
    t0 = time.time()
    net = real.network
    states, counts, spectra, n_unconv, continuum = an.stored_states(real, rng, n_starts)
    out: Dict[str, object] = {"stable_states": len(states), "unsettled_starts": n_unconv}
    if net["type"] == "rotor":
        E = real.V(np.array(states)) if states else np.array([np.inf])
        wsum = float(np.sum(np.abs(real.graph["wa"]) + np.abs(real.graph["wc"])))
        f = float((E.min() + wsum) / wsum)
        g = states[int(np.argmin(E))]
        k = an.kappa_spectrum(real, g)
        out.update(frustration=f, satisfied=bool(f < 1e-6),
                   flat_directions=int(np.sum(np.abs(k.real) < 1e-6 * max(1.0, float(np.abs(k.real).max())))))
    elif net["type"] == "spin":
        E = real.V(np.array(states))
        g = states[int(np.argmin(E))]
        unsat = sum(1 for i, j, s in net["edges"] if s * g[i] * g[j] < 0)
        out.update(unsatisfied_bonds=int(unsat), satisfied=bool(unsat == 0))
    else:
        out.update(oscillates=bool(len(states) == 0 and n_unconv > 0))
    out["seconds"] = time.time() - t0
    return out


def consistent(pred: Dict[str, object], meas: Dict[str, object], kind: str) -> bool:
    if kind == "rotor":
        return (pred["satisfiable"] == meas["satisfied"]
                and pred["continuous_symmetry"] == (meas["flat_directions"] > 0))
    if kind == "spin":
        return pred["satisfiable"] == meas["satisfied"]
    ok = True
    if not pred["positive_loop"]:
        ok &= meas["stable_states"] <= 1
    if pred["monotone"]:
        ok &= not meas["oscillates"]
    return bool(ok)


# ------------------------------------------------------------------------------------------------ random networks
def _connected_edges(rng, n: int, extra: int) -> List[Tuple[int, int]]:
    """A random spanning tree plus `extra` further distinct edges."""
    order = rng.permutation(n)
    edges = {tuple(sorted((int(order[k]), int(order[rng.integers(0, k)])))) for k in range(1, n)}
    pairs = [p for p in itertools.combinations(range(n), 2) if p not in edges]
    rng.shuffle(pairs)
    edges |= set(pairs[:extra])
    return sorted(edges)


def random_rotor_network(rng, n: int, extra: int, kinds=("reflect", "align", "anti"), weights=(0.6, 0.25, 0.15)):
    pos = rng.uniform(0, 3, (n, 2))
    pairs = _connected_edges(rng, n, extra)
    ks = rng.choice(list(kinds), size=len(pairs), p=list(weights))
    return rotor_network(pos, [(i, j, str(k)) for (i, j), k in zip(pairs, ks)])


def random_spin_network(rng, n: int, extra: int, p_anti: float = 0.4):
    pairs = _connected_edges(rng, n, extra)
    return spin_network(n, [(i, j, -1 if rng.random() < p_anti else 1) for i, j in pairs])


def random_gene_network(rng, n: int, p_repress: float = 0.5, max_inputs: int = 2, balanced: bool = False):
    """Each gene has one or two regulators. If `balanced`, the signs follow a random gauge (sign of an edge =
    product of the gauges of its ends), so that every undirected cycle is positive: a monotone network."""
    gauge = rng.choice([-1, 1], size=n)
    regs = []
    for j in range(n):
        k = int(rng.integers(1, max_inputs + 1))
        srcs = rng.choice([i for i in range(n) if i != j], size=k, replace=False)
        for i in srcs:
            s = int(gauge[i] * gauge[j]) if balanced else (-1 if rng.random() < p_repress else 1)
            regs.append((int(i), j, s))
    return gene_network(n, regs)


def benchmark(rng, counts: Dict[str, int], sizes=(4, 5, 6, 7, 8)) -> List[Dict[str, object]]:
    rows = []
    for kind, number in counts.items():
        for t in range(number):
            n = int(sizes[t % len(sizes)])
            if kind == "rotor":
                real = random_rotor_network(rng, n, extra=int(rng.integers(1, n)))
            elif kind == "spin":
                real = random_spin_network(rng, n, extra=int(rng.integers(1, n)))
            else:
                real = random_gene_network(rng, min(n, 6), balanced=(t % 3 == 0))
            pred = predict(real)
            meas = measure(real, rng)
            rows.append({"kind": kind, "n": real.carrier.dim, "edges": real.network["edges"], "prediction": pred,
                         "measured": meas, "consistent": consistent(pred, meas, kind)})
    return rows
