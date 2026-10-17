"""The order of writes: the record measured against the uncoupled reference (the holonomy rule as a measurement).

A structural rule of the memory module says that two writes leave a record of their order only where they fail to
commute on the coordinates that retain the state; the failure is the holonomy of the cycle A B A^-1 B^-1, whose
leading term is the commutator of the two operations. Non-commutation is necessary for a record and is not
sufficient: a device that overwrites its state has writes that do not commute and records only the last of them.
This module turns the rule into a measurement in the same way as the return-point rule is turned into one by
`counterpart`: an index, a unique reference and an excess.

The material is a network of n bistable coordinates q_i in the double well V = (q^2 - w^2)^2 / 4, coupled through
the Laplacian of a graph with strength c, with a response time tau_i per coordinate assigned in breadth-first order
from node 0 (the convention of the composition census of 2026-10-10). A write is a force pulse on a support of
nodes, with a sign, an amplitude, a frequency address nu that sets its loss factor, and a spatial pattern; after
the writes the network relaxes with no force, and the final state is certified as a stable fixed point. All orders
of the named writes are applied; the final states are compared by their sign patterns.

    distinct D       the number of distinct final states over all orders
    split            the number of last writes after which the final state still depends on the earlier order
    uncoupled        the same elements and writes with the coupling removed: each element records the last write
                     that moved it, so this reference has every memory except the one created by the couplings
    gain             log2 D - log2 D_ref, the record of order created by the couplings
    index H          the smallest number of edits, writes removed from a node or couplings cut, after which every
                     write acts on its own independent component and all writes commute exactly (exact up to 10
                     nodes); H = 0 means the writes commute and no order can be recorded

The composition census over every connected graph with 3 to 7 nodes gives the typical record for comparison: with
attractive couplings the share of graphs that record the order falls from 91 % (trees) to 16 % (dense graphs); with
antiferromagnetic couplings it is 26 to 27 % for trees, sparse and medium graphs and 63 % for dense ones.
"""
from __future__ import annotations

import itertools
from collections import deque
from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .spec import SCHEMA, SpecError

KIND = "order"
MAX_UNITS = 12
EXACT_INDEX_UP_TO = 10
DEFAULT_NU = {"A": 0.02, "B": 1.5, "C": float(np.sqrt(0.02 * 1.5))}
DEFAULT_SIGN = {"A": 1.0, "B": -1.0, "C": 1.0}
PULSE = 150 * 0.025
RELEASE = 20000 * 0.025
CENSUS = {  # share of graphs recording the order, by density class, from the census (support N-1, unit amplitude)
    "attractive": {"tree": 0.913, "sparse": 0.789, "medium": 0.456, "dense": 0.156},
    "antiferromagnetic": {"tree": 0.261, "sparse": 0.273, "medium": 0.269, "dense": 0.625},
}


@dataclass
class Write:
    name: str
    support: np.ndarray          # 0/1 mask over the nodes in breadth-first order
    amplitude: float
    sign: float
    nu: float


@dataclass
class OrderNetwork:
    name: str
    n: int
    edges: List[Tuple[int, int]]          # in breadth-first order
    signs: List[float]                    # +1 attractive, -1 antiferromagnetic, per edge
    coupling: float
    well: float
    writes: List[Write]
    pulse: float = PULSE
    release: float = RELEASE
    order: Optional[List[int]] = None     # original node of each breadth-first index

    def __post_init__(self):
        self.tau = np.geomspace(0.08, 40.0, self.n)
        self.drag = 0.25 + self.tau ** 0.35
        A = np.zeros((self.n, self.n))
        for (i, j), s in zip(self.edges, self.signs):
            A[i, j] += s
            A[j, i] += s
        self.L = np.diag(np.abs(A).sum(axis=1)) - A

    # dynamics
    def loss(self, nu: float) -> np.ndarray:
        z = 2 * np.pi * nu * self.tau
        return z / (1 + z * z)

    def force(self, w: Write) -> np.ndarray:
        pattern = 1.0 + 0.12 * np.sin(1.3 * np.arange(self.n) + np.log10(w.nu))
        return w.sign * w.amplitude * self.loss(w.nu) * pattern * w.support

    def gradient(self, q: np.ndarray) -> np.ndarray:
        return q * (q * q - self.well ** 2) + self.coupling * (self.L @ q)

    def hessian(self, q: np.ndarray) -> np.ndarray:
        return np.diag(3 * q * q - self.well ** 2) + self.coupling * self.L

    def advance(self, q: np.ndarray, force: np.ndarray, duration: float) -> np.ndarray:
        from scipy.integrate import solve_ivp
        sol = solve_ivp(lambda t, y: (force - self.gradient(y)) / self.drag, (0.0, duration), q, method="LSODA",
                        rtol=1e-10, atol=1e-12)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y[:, -1]

    def certify(self, q: np.ndarray, max_time: float = 40000.0) -> Optional[np.ndarray]:
        """The stable fixed point the state relaxes to, or None if it is not reached within max_time."""
        from scipy.optimize import root
        elapsed = 0.0
        zero = np.zeros(self.n)
        while elapsed <= max_time:
            if np.max(np.abs(self.gradient(q))) < 1e-7:
                sol = root(lambda x: self.gradient(x), q, jac=lambda x: self.hessian(x), tol=1e-13)
                eq = sol.x
                eig = float(np.linalg.eigvalsh(self.hessian(eq))[0])
                if np.max(np.abs(self.gradient(eq))) < 1e-11 and eig > 1e-8 and np.max(np.abs(eq - q)) < 1e-4:
                    return eq
            q = self.advance(q, zero, 500.0)
            elapsed += 500.0
        return None

    def retained(self, names: Sequence[str]) -> Optional[np.ndarray]:
        """The certified state after the writes in the given order and a field-free release."""
        by = {w.name: w for w in self.writes}
        q = np.full(self.n, -self.well)
        for name in names:
            q = self.advance(q, self.force(by[name]), self.pulse)
        q = self.advance(q, np.zeros(self.n), self.release)
        return self.certify(q)


# ------------------------------------------------------------------------------------------------ specification
def bfs_order(n: int, edges: Sequence[Tuple[int, int]]) -> List[int]:
    adj = [[] for _ in range(n)]
    for i, j in edges:
        adj[i].append(j)
        adj[j].append(i)
    seen, order, dq = {0}, [0], deque([0])
    while dq:
        u = dq.popleft()
        for v in sorted(adj[u]):
            if v not in seen:
                seen.add(v); order.append(v); dq.append(v)
    for v in range(n):
        if v not in seen:
            order.append(v)
    return order


def graph_edges(spec_graph, n: int) -> List[Tuple[int, int]]:
    if isinstance(spec_graph, str):
        if spec_graph == "ring":
            return [(i, (i + 1) % n) for i in range(n)]
        if spec_graph == "chain":
            return [(i, i + 1) for i in range(n - 1)]
        if spec_graph == "complete":
            return [(i, j) for i in range(n) for j in range(i + 1, n)]
        if spec_graph == "star":
            return [(0, j) for j in range(1, n)]
        raise SpecError(f"graph must be ring, chain, complete, star or a list of edges, not {spec_graph!r}")
    edges = []
    for e in spec_graph:
        i, j = int(e[0]), int(e[1])
        if not (0 <= i < n and 0 <= j < n) or i == j:
            raise SpecError(f"edge {e}: nodes are numbered 0 to {n - 1} and an edge joins two different nodes")
        edges.append((i, j))
    return edges


def from_spec(spec: Dict) -> OrderNetwork:
    """The network of a specification of kind ``order``; SpecError says what to change."""
    if spec.get("schema") != SCHEMA:
        raise SpecError(f"schema must be {SCHEMA!r}")
    if spec.get("kind") != KIND:
        raise SpecError(f"kind must be {KIND!r} for this command")
    o = spec.get("order")
    if not isinstance(o, dict):
        raise SpecError("A specification of kind 'order' needs an 'order' object")
    n = int(o.get("units", 0))
    if not 2 <= n <= MAX_UNITS:
        raise SpecError(f"units must be between 2 and {MAX_UNITS}")
    edges = graph_edges(o.get("graph", "ring"), n)
    sign = o.get("couplings", "attractive")
    if sign not in ("attractive", "antiferromagnetic"):
        raise SpecError("couplings must be 'attractive' or 'antiferromagnetic'")
    ordered = o.get("graph") == "ring"            # the ring keeps its own order, as the census control did
    perm = list(range(n)) if ordered else bfs_order(n, edges)
    pos = {v: k for k, v in enumerate(perm)}
    bedges = [(pos[i], pos[j]) for i, j in edges]
    writes_spec = o.get("writes")
    if not isinstance(writes_spec, dict) or not 2 <= len(writes_spec) <= 4:
        raise SpecError("writes must be an object of two to four named writes, e.g. {'A': {...}, 'B': {...}}")
    writes = []
    for name, w in writes_spec.items():
        w = w or {}
        mask = np.zeros(n)
        if "support" in w:
            for v in w["support"]:
                if not 0 <= int(v) < n:
                    raise SpecError(f"write {name}: support nodes are numbered 0 to {n - 1}")
                mask[pos[int(v)]] = 1.0
        else:
            k = int(w.get("k", n - 1))
            if not 1 <= k <= n:
                raise SpecError(f"write {name}: k must be between 1 and {n}")
            lo = (n - k) // 2
            mask[lo:lo + k] = 1.0
        writes.append(Write(name=str(name), support=mask, amplitude=float(w.get("amplitude", 1.0)),
                            sign=float(w.get("sign", DEFAULT_SIGN.get(name, 1.0))),
                            nu=float(w.get("frequency", DEFAULT_NU.get(name, 0.1)))))
    return OrderNetwork(name=str(spec.get("name", "order")), n=n, edges=bedges,
                        signs=[1.0 if sign == "attractive" else -1.0] * len(bedges),
                        coupling=float(o.get("coupling", 0.13)), well=float(o.get("well", 0.65)), writes=writes,
                        pulse=float(o.get("pulse", PULSE)), release=float(o.get("release", RELEASE)), order=perm)


# ------------------------------------------------------------------------------------------------ measurement
def _signs(q: np.ndarray) -> str:
    return "".join("+" if v > 0 else "-" for v in q)


def measure(net: OrderNetwork) -> Dict[str, object]:
    """Final states for every order of the writes, their number, and the split: last writes after which the state
    still depends on the earlier order."""
    names = [w.name for w in net.writes]
    states: Dict[str, Optional[str]] = {}
    for perm in itertools.permutations(names):
        eq = net.retained(perm)
        states["".join(perm)] = None if eq is None else _signs(eq)
    resolved = {o: s for o, s in states.items() if s is not None}
    by_last: Dict[str, set] = {}
    for o, s in resolved.items():
        by_last.setdefault(o[-1], set()).add(s)
    distinct = len(set(resolved.values()))
    return {"states": states, "unresolved": len(states) - len(resolved), "distinct": distinct,
            "capacity": float(np.log2(max(distinct, 1))),
            "split": sum(1 for v in by_last.values() if len(v) > 1), "last_writes": len(by_last)}


def uncoupled(net: OrderNetwork) -> OrderNetwork:
    """The same elements and writes with the coupling removed: the last-write reference."""
    return replace(net, name=f"{net.name} (uncoupled reference)", coupling=0.0)


def commutation_index(net: OrderNetwork, exact_up_to: int = EXACT_INDEX_UP_TO) -> Dict[str, object]:
    """The smallest number of edits after which the writes commute exactly: a write removed from a node, or a
    coupling cut, such that every connected component carries the writes of one name only. Exact by enumeration of
    the assignments of nodes to writes up to exact_up_to nodes; above, a greedy assignment, an upper bound."""
    names = [w.name for w in net.writes]
    on = [{w.name for w in net.writes if w.support[v] > 0} for v in range(net.n)]
    m = len(names)

    def cost(assign):
        c = sum(len(on[v]) - (1 if names[assign[v]] in on[v] else 0) for v in range(net.n))
        c += sum(1 for i, j in net.edges if assign[i] != assign[j])
        return c
    if net.n <= exact_up_to:
        best, best_assign = None, None
        for assign in itertools.product(range(m), repeat=net.n):
            c = cost(assign)
            if best is None or c < best:
                best, best_assign = c, assign
        exact = True
    else:
        assign = [max(range(m), key=lambda k: names[k] in on[v]) for v in range(net.n)]
        improved = True
        while improved:
            improved = False
            for v in range(net.n):
                for k in range(m):
                    trial = list(assign); trial[v] = k
                    if cost(trial) < cost(assign):
                        assign, improved = trial, True
        best, best_assign, exact = cost(assign), tuple(assign), False
    return {"H": int(best), "exact": exact, "assignment": {names[k]: [v for v in range(net.n) if best_assign[v] == k] for k in range(m)}}


def density_class(net: OrderNetwork) -> str:
    e, n = len(net.edges), net.n
    if e == n - 1:
        return "tree"
    m = (e - (n - 1)) / max(1, n * (n - 1) // 2 - (n - 1))
    return "sparse" if m < 0.34 else "medium" if m < 0.67 else "dense"


def construct(net: OrderNetwork, exact_up_to: int = EXACT_INDEX_UP_TO) -> Dict[str, object]:
    """The order record of the network, its uncoupled reference, the gain, the commutation index and the typical
    record of the census for the network's density and coupling sign."""
    idx = commutation_index(net, exact_up_to)
    meas = measure(net)
    ref = measure(uncoupled(net))
    kind = "attractive" if all(s > 0 for s in net.signs) else "antiferromagnetic" if all(s < 0 for s in net.signs) else "mixed"
    cls = density_class(net)
    out: Dict[str, object] = {
        "rule": "the order of writes is recorded only where the writes fail to commute on the retaining coordinates "
                "(necessary, not sufficient)",
        "index": idx, "network": meas, "uncoupled_reference": ref,
        "gain": meas["capacity"] - ref["capacity"], "split_gain": meas["split"] - ref["split"],
        "density_class": cls, "couplings": kind,
        "typical_record_in_census": CENSUS.get(kind, {}).get(cls)}
    out["statement"] = statement(out)
    return out


def statement(r: Dict[str, object]) -> str:
    idx, net, ref = r["index"], r["network"], r["uncoupled_reference"]
    parts = []
    if idx["H"] == 0:
        parts.append("The writes commute (H = 0): no order can be recorded.")
    else:
        parts.append(f"The writes fail to commute: {idx['H']} edits ({'exact' if idx['exact'] else 'bound'}) "
                     "would make them commute.")
    parts.append(f"Over the {len(net['states'])} orders the network ends in {net['distinct']} distinct states "
                 f"({net['capacity']:.2f} bits) and the order after the same last write is retained in "
                 f"{net['split']} of {net['last_writes']} cases; the uncoupled reference ends in {ref['distinct']} "
                 f"({ref['capacity']:.2f} bits, split {ref['split']}): the couplings add {r['gain']:.2f} bits.")
    typ = r.get("typical_record_in_census")
    if typ is not None:
        parts.append(f"In the census, {100 * typ:.0f} % of {r['density_class']} graphs with {r['couplings']} "
                     "couplings record the order.")
    if net["unresolved"]:
        parts.append(f"{net['unresolved']} orders did not reach a certified state.")
    return " ".join(parts)
