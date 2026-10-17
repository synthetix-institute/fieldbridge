"""Counterpart constructions: the smallest change that restores a structural rule, and the excess measured against it.

A structural rule says that a property of a material follows from the signs or the algebra of its equations. When a
material violates the rule, the rule alone says only that the property is no longer guaranteed. The counterpart
turns this into a measurement. It is the material obtained by the smallest change of one layer of the description
that makes the rule hold, with everything else kept: the same graph, magnitudes, fields and thresholds. The
difference between a quantity measured on the material and on its counterpart is the consequence of the violation,
and it is a number an experiment can compare.

The first rule treated here is the return to a turning point of interacting hysterons (fieldbridge.memory.hysterons).
Count the drive as one more element, joined to every driven element with the sign of its drive. If every loop of this
signed network can be satisfied at once (the network is balanced in the sense of Harary, Michigan Math. J. 2, 143,
1953), a relabeling makes every coupling cooperative and the drive uniform, and the return to every turning point is
exact. The frustration index F is the smallest number of edges, couplings or actions of the drive, whose signs must be
reversed to make the network balanced; a relabeling g that attains it leaves exactly F unsatisfiable edges. Reversing
these edges gives the balanced counterpart: same couplings in magnitude, same fields, thresholds and strengths of the
drive, and return-point memory by construction.

Measured against the counterpart, the number of distinguishable states of the network separates into two parts:

    loss = C_bal - C_0     the counterpart's locking below its uncoupled elements (C_0: all couplings removed)
    gain = C - C_bal       states that exist only because of the unsatisfiable edges

with C = log2 of the number of states that some history of the drive reaches from negative saturation. States that
need three or more reversals of the drive are reported as deep memory. The return test is exhaustive: every reachable
state is used as a turning point, with every excursion that does not pass the previous turning point.

The transition graph is enumerated exactly (avalanches relaxed least stable element first, as in Run); its size grows
quickly with the number of elements, so the enumeration stops at max_states and the result says so. The frustration
index is exact up to exact_up_to elements (2^n relabelings) and above that is the best of local searches, an upper
bound.
"""
from __future__ import annotations

import itertools
from collections import deque
from dataclasses import replace
from typing import Dict, List, Optional, Tuple

import numpy as np

from .hysterons import AvalancheError, Hysterons, predict

EXACT_UP_TO = 16          # elements up to which the frustration index is computed over all relabelings
MAX_STATES = 20000        # states of a transition graph before the enumeration stops
EPS = 1e-9                # the drive is moved this far past an event, in units of the largest field
DEEP = 3                  # reversals of the drive from which a state counts as deep memory


class TooManyStates(RuntimeError):
    """The transition graph has more states than the enumeration is allowed to visit."""


# ------------------------------------------------------------------------------------------------ signed network
def _edges(model: Hysterons) -> Tuple[List[Tuple[int, int, float]], List[Tuple[int, float]]]:
    """The edges of the network with the drive: couplings (i, j, sign), each reciprocal coupling once and a coupling
    that is not reciprocal once per direction, and actions of the drive (i, sign of eta_i)."""
    couplings = [(i, j, float(np.sign(J))) for i, j, J in model.undirected() if J != 0]
    drive = [(int(i), float(np.sign(model.eta[i]))) for i in np.flatnonzero(model.eta != 0)]
    return couplings, drive


def _unsatisfied(g: np.ndarray, couplings, drive) -> Tuple[List[int], List[int]]:
    bad_c = [k for k, (i, j, s) in enumerate(couplings) if g[i] * s * g[j] < 0]
    bad_d = [k for k, (i, s) in enumerate(drive) if g[i] * s < 0]
    return bad_c, bad_d


def frustration(model: Hysterons, rng: Optional[np.random.Generator] = None, exact_up_to: int = EXACT_UP_TO,
                restarts: int = 40) -> Dict[str, object]:
    """Frustration index of the signed network of the elements and the drive (the drive has the fixed label +1).

    Returns F, a relabeling g that attains it (or the best one found), whether F is exact, the unsatisfiable couplings
    and actions of the drive under g, and the numbers of edges. Up to exact_up_to elements every relabeling is
    tried; above, F is the smallest count of `restarts` greedy local searches started from random relabelings, an
    upper bound on the index."""
    couplings, drive = _edges(model)
    n = model.n
    if n <= exact_up_to:
        G = np.array(list(itertools.product([1.0, -1.0], repeat=n)))
        count = np.zeros(len(G), dtype=int)
        if couplings:
            ci = np.array([c[0] for c in couplings]); cj = np.array([c[1] for c in couplings])
            cs = np.array([c[2] for c in couplings])
            count += (G[:, ci] * G[:, cj] * cs < 0).sum(axis=1)
        if drive:
            di = np.array([d[0] for d in drive]); ds = np.array([d[1] for d in drive])
            count += (G[:, di] * ds < 0).sum(axis=1)
        best = int(np.argmin(count))
        g, exact = G[best], True
    else:
        rng = rng if rng is not None else np.random.default_rng(0)
        incident: List[List[Tuple[int, float]]] = [[] for _ in range(n)]
        for i, j, s in couplings:
            incident[i].append((j, s))
            incident[j].append((i, s))
        field = np.zeros(n)
        for i, s in drive:
            field[i] += s
        g, best_count = None, None
        for _ in range(restarts):
            trial = rng.choice([-1.0, 1.0], n)
            improved = True
            while improved:
                improved = False
                for i in range(n):
                    # satisfied minus unsatisfied edges at i; flipping g_i reverses its sign
                    balance = sum(s * trial[i] * trial[j] for j, s in incident[i]) + field[i] * trial[i]
                    if balance < 0:
                        trial[i] = -trial[i]
                        improved = True
            bad_c, bad_d = _unsatisfied(trial, couplings, drive)
            if best_count is None or len(bad_c) + len(bad_d) < best_count:
                g, best_count = trial.copy(), len(bad_c) + len(bad_d)
        exact = False
    bad_c, bad_d = _unsatisfied(g, couplings, drive)
    return {"F": len(bad_c) + len(bad_d), "exact": exact, "relabeling": [int(x) for x in g],
            "unsatisfied_couplings": [(int(couplings[k][0]), int(couplings[k][1])) for k in bad_c],
            "unsatisfied_drive": [int(drive[k][0]) for k in bad_d],
            "edges_couplings": len(couplings), "edges_drive": len(drive)}


def balanced_counterpart(model: Hysterons, relabeling) -> Hysterons:
    """The material with the edges reversed that are unsatisfied under the relabeling: couplings g_i J_ij g_j < 0
    (both directions of a reciprocal coupling) and actions of the drive g_i eta_i < 0. Magnitudes, fields, half-widths
    and the graph are unchanged."""
    g = np.asarray(relabeling, dtype=float)
    bonds = [(i, j, -J if g[i] * J * g[j] < 0 else J) for i, j, J in model.bonds]
    eta = np.where(g * model.eta < 0, -model.eta, model.eta)
    return replace(model, name=f"{model.name} (balanced counterpart)", bonds=bonds, eta=eta)


def uncoupled(model: Hysterons) -> Hysterons:
    """The same elements, fields, half-widths and drive with every coupling removed (independent hysterons)."""
    return replace(model, name=f"{model.name} (uncoupled)", bonds=[])


# ------------------------------------------------------------------------------------------------ transition graph
def _matrix(model: Hysterons) -> np.ndarray:
    J = np.zeros((model.n, model.n))
    for i, j, w in model.bonds:
        J[i, j] += w
    return J


class _Graph:
    """Stable states of a material and the transitions between them under a slow drive.

    A stable state is stable on an interval of the drive; the transition up leads to the state after the avalanche
    just above the upper end of that interval, the transition down to the state after the avalanche just below the
    lower end. Avalanches are relaxed least stable element first."""

    def __init__(self, model: Hysterons):
        self.m = model
        self.J = _matrix(model)
        self.eps = EPS * float(np.abs(self.J).sum(axis=1).max(initial=0.0) + np.abs(model.h).max(initial=0.0)
                               + model.b.max(initial=0.0) + 1.0)

    def field(self, s: np.ndarray, H: float) -> np.ndarray:
        return self.J @ s + self.m.h + self.m.eta * H

    def interval(self, s: np.ndarray) -> Optional[Tuple[float, float]]:
        a = s * self.field(s, 0.0) + self.m.b
        c = s * self.m.eta
        if np.any(a[c == 0] < 0):
            return None
        lo = float(np.max(-a[c > 0] / c[c > 0])) if np.any(c > 0) else -np.inf
        hi = float(np.min(a[c < 0] / -c[c < 0])) if np.any(c < 0) else np.inf
        return (lo, hi) if lo <= hi else None

    def relax(self, s: np.ndarray, H: float) -> np.ndarray:
        s = s.copy()
        f = self.field(s, H)
        for _ in range(50 * self.m.n + 1000):
            st = s * f + self.m.b
            i = int(np.argmin(st))
            if st[i] >= 0:
                return s
            s[i] = -s[i]
            f += 2.0 * self.J[:, i] * s[i]
        raise AvalancheError("an avalanche did not end; with couplings that are not reciprocal the dynamics can cycle")

    def step(self, s: np.ndarray, up: bool) -> Optional[Tuple[np.ndarray, float]]:
        lo, hi = self.interval(s)
        edge = hi if up else lo
        if not np.isfinite(edge):
            return None
        H = edge + self.eps if up else edge - self.eps
        return self.relax(s, H), H

    def saturated(self) -> np.ndarray:
        s = np.where(self.m.eta > 0, -1.0, np.where(self.m.eta < 0, 1.0, -1.0))
        return self.relax(s, self.m.saturation())


def _key(s: np.ndarray) -> bytes:
    return (s > 0).astype(np.uint8).tobytes()


def reachable(model: Hysterons, max_states: int = MAX_STATES) -> Dict[str, object]:
    """The states reached from negative saturation by any history of the drive, with their up and down successors.
    Raises TooManyStates beyond max_states and AvalancheError if an avalanche does not end."""
    graph = _Graph(model)
    s0 = graph.saturated()
    states, index = [s0], {_key(s0): 0}
    up: List[Optional[int]] = [None]
    down: List[Optional[int]] = [None]
    queue = deque([0])
    while queue:
        k = queue.popleft()
        for direction, edges in ((True, up), (False, down)):
            nxt = graph.step(states[k], direction)
            if nxt is None:
                continue
            t = _key(nxt[0])
            if t not in index:
                if len(states) >= max_states:
                    raise TooManyStates(f"more than {max_states} states")
                index[t] = len(states)
                states.append(nxt[0])
                up.append(None)
                down.append(None)
                queue.append(index[t])
            edges[k] = index[t]
    return {"states": states, "up": up, "down": down, "count": len(states), "graph": graph}


def reversals(tgraph: Dict[str, object]) -> List[int]:
    """For each reachable state, the smallest number of reversals of the drive in a history from negative saturation
    that ends in it: moving on in the current direction costs nothing, a change of direction costs one."""
    n = tgraph["count"]
    inf = 1 << 30
    dist = {(i, d): inf for i in range(n) for d in (True, False)}
    dist[(0, True)] = 0
    dq = deque([(0, True)])
    while dq:
        i, d = dq.popleft()
        c = dist[(i, d)]
        nxt = (tgraph["up"] if d else tgraph["down"])[i]
        if nxt is not None and dist[(nxt, d)] > c:
            dist[(nxt, d)] = c
            dq.appendleft((nxt, d))
        if dist[(i, not d)] > c + 1:
            dist[(i, not d)] = c + 1
            dq.append((i, not d))
    return [min(dist[(i, True)], dist[(i, False)]) for i in range(n)]


def exact_return(tgraph: Dict[str, object]) -> Dict[str, int]:
    """Return to a turning point, tested on every reachable state.

    A history that turns upward in a state M (reached by a falling drive, or saturation) rises along the up path of
    M. Every state A on that path can hold the next turning point H1, anywhere between the field at which A was
    created and the upper end of its interval. From A the drive goes down by any number of avalanches without passing
    the previous turning point and comes back up to H1; the state reached must be A for every admissible H1. The same
    with the directions exchanged. Counts the excursions tested and those that fail for some H1."""
    graph: _Graph = tgraph["graph"]
    states = tgraph["states"]
    tested = failed = 0
    starts = {True: {0} | {t for t in tgraph["down"] if t is not None},
              False: {t for t in tgraph["up"] if t is not None}}
    for rising, members in starts.items():
        r = 1.0 if rising else -1.0

        def ends(S):
            lo, hi = graph.interval(S)
            return (lo, hi) if rising else (-hi, -lo)

        for mi in sorted(members):
            M = states[mi]
            loM, hiM = ends(M)
            bound = loM + min(graph.eps, 0.5 * (hiM - loM)) if np.isfinite(loM) else -np.inf
            A, created = M, bound
            for _ in range(graph.m.n + 2):
                lo, hi = ends(A)
                if not np.isfinite(hi):
                    break
                B = A
                while True:
                    d = graph.step(B, not rising)
                    if d is None or r * d[1] <= bound:
                        break
                    B = d[0]
                    C, prev_hi, bad = B, -np.inf, False
                    while True:
                        _, hiC = ends(C)
                        if hiC >= created and prev_hi < hi and not np.array_equal(C, A):
                            bad = True
                            break
                        if hiC >= hi:
                            break
                        prev_hi = hiC
                        C = graph.step(C, rising)[0]
                    tested += 1
                    failed += int(bad)
                nxt = graph.step(A, rising)
                if nxt is None:
                    break
                A, created = nxt[0], r * nxt[1]
    return {"tested": tested, "failed": failed}


def measure(model: Hysterons, max_states: int = MAX_STATES, return_test: bool = True) -> Dict[str, object]:
    """Distinguishable states of one material: count, C = log2(count), deep memory (states that need DEEP or more
    reversals of the drive, in bits above the shallow count), the largest number of reversals needed, and the
    exhaustive return test. An enumeration that stops reports `censored`; an avalanche that does not end, `cycle`."""
    try:
        tg = reachable(model, max_states)
    except TooManyStates:
        return {"censored": True, "max_states": max_states}
    except AvalancheError:
        return {"cycle": True}
    depth = np.array(reversals(tg))
    shallow = int((depth < DEEP).sum())
    out: Dict[str, object] = {"count": int(tg["count"]), "C": float(np.log2(tg["count"])),
                              "deep": float(np.log2(tg["count"]) - np.log2(max(shallow, 1))),
                              "max_reversals": int(depth.max())}
    if return_test:
        out["return"] = exact_return(tg)
    return out


# ------------------------------------------------------------------------------------------------ the construction
def construct(model: Hysterons, rng: Optional[np.random.Generator] = None, max_states: int = MAX_STATES,
              exact_up_to: int = EXACT_UP_TO, return_test: bool = True) -> Dict[str, object]:
    """The balanced counterpart of a hysteron material and the measurement of the material against it."""
    fr = frustration(model, rng, exact_up_to=exact_up_to)
    cp = balanced_counterpart(model, fr["relabeling"])
    check = predict(cp)
    if not check["guaranteed"]:
        raise RuntimeError("the counterpart is not balanced: reversing the unsatisfied edges must balance the network")
    net, bal = measure(model, max_states, return_test), measure(cp, max_states, return_test)
    free = measure(uncoupled(model), max_states, return_test=False)
    out: Dict[str, object] = {
        "rule": "return to a turning point (R8): guaranteed when the network of couplings with the drive is balanced",
        "frustration": fr, "counterpart": {"name": cp.name, "reversed_couplings": fr["unsatisfied_couplings"],
                                           "reversed_drive": fr["unsatisfied_drive"]},
        "network": net, "balanced": bal, "uncoupled": free}
    if all("C" in x for x in (net, bal, free)):
        out["gain"] = net["C"] - bal["C"]
        out["loss"] = bal["C"] - free["C"]
        out["deep_gain"] = net["deep"] - bal["deep"]
    out["statement"] = statement(out)
    return out


def statement(result: Dict[str, object]) -> str:
    """The result in words."""
    fr = result["frustration"]
    F = fr["F"]
    index = f"F = {F}" if fr["exact"] else f"F ≤ {F} (local search)"
    if F == 0:
        return ("The network with the drive is balanced (F = 0): the material is its own counterpart, and its return "
                "to every turning point is guaranteed.")
    parts = [f"Reversing {F} unsatisfiable edges ({len(fr['unsatisfied_couplings'])} couplings, "
             f"{len(fr['unsatisfied_drive'])} actions of the drive; {index}) gives the balanced counterpart, whose "
             "return to every turning point is guaranteed."]
    net, bal = result["network"], result["balanced"]
    if "gain" in result:
        parts.append(f"The material reaches {net['count']} states ({net['C']:.2f} bits), its counterpart "
                     f"{bal['count']} ({bal['C']:.2f} bits): a gain of {result['gain']:.2f} bits from the "
                     f"unsatisfiable edges, against a loss of {-result['loss']:.2f} bits of the counterpart below its "
                     "uncoupled elements.")
    else:
        parts.append("The enumeration of states stopped before the end; the gain is not measured.")
    if "return" in net and "return" in bal:
        parts.append(f"Return test: the material fails {net['return']['failed']} of {net['return']['tested']} "
                     f"excursions, the counterpart {bal['return']['failed']} of {bal['return']['tested']}.")
    return " ".join(parts)
