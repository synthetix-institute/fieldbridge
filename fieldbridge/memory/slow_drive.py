"""Return to a turning point for a material given by equations: its transition graph under a slow drive, the exact
return test on that graph, and the reduction to interacting hysterons (specification kind ``equations`` with a block
``slow_drive``).

The hysteron description (fieldbridge.memory.hysterons) starts from elements that switch at thresholds. A material
given by equations has no thresholds written in it: it has a drift F(q; H) and a parameter H through which a slow
drive acts. Under a quasi-static drive and without noise its state sits on a stable equilibrium that moves with H
until the equilibrium disappears at a fold; there the state relaxes along the drift to another stable equilibrium.
This is the dynamics that the hysteron model assumes, and it is followed here on the equations themselves:

  transition graph   every stable equilibrium reached from the low end of the range of H by some history of the
                     drive, the interval of H on which it exists (between its folds, or up to an end of the range),
                     and the equilibrium reached by relaxation just past each fold. The number of these states, the
                     reversals of the drive that each needs and the exhaustive return test of
                     fieldbridge.memory.counterpart are computed on this graph. They are therefore statements about
                     the equations under a quasi-static drive, with no model of elements in between.
  elements           observables o_i(q) whose signs label the states, for coupled bistable units a coordinate of
                     each unit. They serve the reduction only: a state is identified by its equilibrium, and two
                     equilibria with the same signs are two states.
  hysteron model     at a fold, the element whose observable changes sign first during the relaxation is the one
                     that lost stability. Each such event is one linear condition
                         s_i (sum_j J_ij s_j + h_i + eta_i H) + b_i = 0
                     on the fields h_i, half-widths b_i and couplings J_ij, where eta_i = +1 or -1 is the direction in
                     which a rising drive switches element i. The conditions are solved by least squares. The residual
                     says how far the elements interact pairwise; the signs of J and eta give the network to which the
                     structural rule of the return applies (hysterons.predict: if no loop through the drive is
                     frustrated, every return is exact). The reduced model is measured too, and its count of states
                     and its return are compared with those of the equations.

Limits. The range of H must saturate the material: one stable equilibrium at each end. A state that does not settle
after a fold (an oscillation, an escape) stops the construction with that reason. So does a degenerate fold, at which
two elements lose stability together and the state reached depends on the direction of an infinitesimal
perturbation: identical elements do this, and a spread of their parameters lifts it. An element whose observable
changes sign along an equilibrium, without a fold, switches continuously; it is reported and is not a hysteron. All
statements are about the noise-free, quasi-static dynamics. Noise lets the state leave an equilibrium before its
fold, which narrows the loops and can change the transition graph.
"""
from __future__ import annotations

import bisect
from collections import deque
from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import analysis as an
from .identity import Realization
from .spec import SpecError, parse_expression, _sympy

STEPS = 400            # continuation steps across the range of the drive
FOLD_TOL = 1e-7        # a fold is located to this fraction of the range
PAST = 1e-4            # after a fold the drive is moved this fraction of the range past it
PERTURB = 1e-5         # size of the perturbations that test a fold for degeneracy, in units of the carrier scale
MAX_TIME = 20000.0     # longest relaxation after a fold, in units of time of the drift
MAX_STATES = 2000      # states of a transition graph before the enumeration stops


class ReductionError(RuntimeError):
    """The quasi-static dynamics is not that of hysterons; the message says why."""


class TooManyStates(RuntimeError):
    pass


@dataclass
class Branch:
    """A stable equilibrium followed along the drive: samples (H increasing), its ends and the signs of the
    elements."""
    H: List[float]
    q: List[np.ndarray]
    lo: float
    hi: float
    signs: np.ndarray
    continuous: List[int] = field(default_factory=list)


def compile_observables(real: Realization, expressions: Sequence[str]):
    """o(q) -> the values of the element observables, expressions in the variables and parameters of the spec."""
    sp = _sympy()
    spec = real.spec
    variables = list(spec["carrier"]["variables"])
    params = list(spec["parameters"])
    vs = [sp.Symbol(v, real=True) for v in variables]
    ps = [sp.Symbol(p, real=True) for p in params]
    names = {**dict(zip(variables, vs)), **dict(zip(params, ps))}
    parsed = [parse_expression(e, names) for e in expressions]
    fn = sp.lambdify(vs + ps, parsed, modules="numpy")

    def obs(q: np.ndarray, p: Dict[str, float]) -> np.ndarray:
        return np.array([float(v) for v in fn(*np.asarray(q, float), *[p[k] for k in params])])

    return obs


def from_spec(spec: Dict, real: Realization) -> Tuple[str, Tuple[float, float], List[str]]:
    """The drive parameter, its range and the element observables of a block ``slow_drive``."""
    block = spec.get("slow_drive")
    if not isinstance(block, dict):
        raise SpecError("A material given by equations needs a block 'slow_drive' with 'parameter', 'range' and "
                        "'elements' for the return to a turning point")
    name = block.get("parameter")
    if name not in real.params:
        raise SpecError(f"slow_drive.parameter {name!r} must be a declared parameter")
    rng_ = block.get("range")
    if not (isinstance(rng_, list) and len(rng_) == 2 and all(isinstance(v, (int, float)) for v in rng_)
            and rng_[0] < rng_[1]):
        raise SpecError("slow_drive.range must be [low, high] with low < high")
    elements = block.get("elements")
    if not (isinstance(elements, list) and elements and all(isinstance(e, str) for e in elements)):
        raise SpecError("slow_drive.elements must be a list of expressions, one per element, whose sign is its state")
    if real.null_modes is not None:
        raise SpecError("a material with a continuous symmetry has no isolated equilibria to follow")
    return str(name), (float(rng_[0]), float(rng_[1])), list(elements)


class EquationGraph:
    """Stable equilibria of the equations and the transitions between them under a slow drive, with the interface
    that counterpart.reversals and counterpart.exact_return use: states are arrays [k] (k the index of a branch),
    interval(S) the ends of its interval, step(S, up) the state after the fold and the drive just past it."""

    def __init__(self, real: Realization, drive: str, H_range: Tuple[float, float], elements: Sequence[str],
                 steps: int = STEPS, seed: int = 0):
        self.real, self.drive = real, drive
        self.lo, self.hi = float(H_range[0]), float(H_range[1])
        span = self.hi - self.lo
        self.dH = span / steps
        self.fold_tol = FOLD_TOL * span
        self.past = PAST * span
        self.eps = self.past
        self.obs = compile_observables(real, elements)
        self.elements = list(elements)
        self.n_elements = len(elements)
        length = an.base_length(real) * np.sqrt(real.carrier.dim)
        self.jump = 0.2 * length
        self.merge = an._merge_tol(real)
        self.perturb = PERTURB * an.base_length(real)
        self.dt = an.safe_dt(real)
        self.rng = np.random.default_rng(seed)
        self.branches: List[Branch] = []
        self.events: List[Dict[str, object]] = []
        self._steps: Dict[Tuple[int, bool], Optional[Tuple[np.ndarray, float]]] = {}
        self.m = SimpleNamespace(n=1)

    # -- equilibria ------------------------------------------------------------------------------------------
    def params(self, H: float) -> Dict[str, float]:
        return {**self.real.params, self.drive: H}

    def stable(self, q: np.ndarray, H: float) -> Optional[np.ndarray]:
        """The equilibrium that Newton's method reaches from q at H, if it is stable."""
        over = {self.drive: H}
        q2, ok = an.newton(self.real, q, over)
        if not ok or not np.all(np.isfinite(q2)):
            return None
        k = an.kappa_spectrum(self.real, q2, over).real
        if k.min() <= 1e-9 * max(1.0, float(np.abs(k).max())):
            return None
        return q2

    def stable_states(self, H: float, n_starts: int = 48) -> List[np.ndarray]:
        """Stable equilibria at H, by relaxation from random starts over three times the carrier scale."""
        starts = self.real.carrier.sample(self.rng, n_starts)
        if self.real.carrier.kind == "euclid":
            starts = 3.0 * starts
        q, _ = an.relax(self.real, starts, {self.drive: H}, max_time=2000.0)
        found: List[np.ndarray] = []
        for row in q:
            s = self.stable(row, H)
            if s is not None and not any(self.real.carrier.distance(s, f) < self.merge for f in found):
                found.append(s)
        return found

    def _follow(self, q: np.ndarray, H: float, direction: int) -> Tuple[List[float], List[np.ndarray], float]:
        """Continue the stable equilibrium from (H, q) until it disappears (the last H at which it exists, located
        to fold_tol) or the range ends (+-inf)."""
        end = self.hi if direction > 0 else self.lo
        Hs, qs, dH = [H], [q], self.dH
        while True:
            if direction * (end - H) <= 0:
                return Hs, qs, direction * np.inf
            H_new = H + direction * min(dH, abs(end - H))
            q_new = self.stable(q, H_new)
            if q_new is not None and self.real.carrier.distance(q_new, q) <= self.jump:
                H, q = H_new, q_new
                Hs.append(H)
                qs.append(q)
                dH = min(2.0 * dH, self.dH)
                continue
            if dH <= self.fold_tol:
                return Hs, qs, H
            dH *= 0.5

    def _branch(self, q: np.ndarray, H: float) -> int:
        H_up, q_up, hi = self._follow(q, H, +1)
        H_dn, q_dn, lo = self._follow(q, H, -1)
        Hs = H_dn[::-1] + H_up[1:]
        qs = q_dn[::-1] + q_up[1:]
        signs = np.array([np.sign(self.obs(x, self.params(h))) for h, x in zip(Hs, qs)])
        middle = signs[len(signs) // 2]
        continuous = [int(i) for i in range(self.n_elements) if len(set(signs[:, i])) > 1]
        self.branches.append(Branch(Hs, qs, lo, hi, middle, continuous))
        return len(self.branches) - 1

    def point(self, k: int, H: float) -> Optional[np.ndarray]:
        """The equilibrium of branch k at H, from its nearest sample."""
        b = self.branches[k]
        if not (b.lo - self.fold_tol <= H <= b.hi + self.fold_tol):
            return None
        i = min(max(bisect.bisect_left(b.H, H), 0), len(b.H) - 1)
        if i > 0 and abs(b.H[i - 1] - H) < abs(b.H[i] - H):
            i -= 1
        return self.stable(b.q[i], H)

    def identify(self, q: np.ndarray, H: float) -> int:
        for k in range(len(self.branches)):
            p = self.point(k, H)
            if p is not None and self.real.carrier.distance(p, q) < self.merge:
                return k
        if len(self.branches) >= MAX_STATES:
            raise TooManyStates(f"more than {MAX_STATES} states")
        return self._branch(q, H)

    # -- the fold ----------------------------------------------------------------------------------------------
    def _settle(self, q_ref: np.ndarray, delta: np.ndarray, H: float) -> Tuple[np.ndarray, List[int]]:
        """Relax without noise from q_ref + delta at H; the stable equilibrium reached and the elements in the order
        in which their observables changed sign. Where the relaxation stops at a saddle (an invariant subspace of
        identical elements does this), it leaves along the unstable direction on the side of the initial
        perturbation, so that opposite perturbations reveal a choice that the drive does not make."""
        over, p = {self.drive: H}, self.params(H)
        o0 = np.sign(self.obs(q_ref, p))
        q, order = self.real.carrier.wrap(q_ref + delta), []
        steps = int(MAX_TIME / self.dt)
        for k in range(steps):
            f = self.real.F(q, **over)
            q = self.real.carrier.wrap(q + self.dt * f)
            o = np.sign(self.obs(q, p))
            for i in np.flatnonzero(o != o0):
                if int(i) not in order:
                    order.append(int(i))
            if k % 20 == 0 and np.linalg.norm(f) < 1e-7:
                s = self.stable(q, H)
                if s is not None:
                    o = np.sign(self.obs(s, p))
                    order += [int(i) for i in np.flatnonzero(o != o0) if int(i) not in order]
                    return s, order
                w, V = np.linalg.eig(an.jacobian(self.real, q, over))
                v = np.real(V[:, int(np.argmax(w.real))])
                side = float(np.dot(v, self.real.carrier.diff(q_ref, q)))
                v *= (1.0 if side >= 0 else -1.0) * 100.0 * self.perturb / max(np.linalg.norm(v), 1e-300)
                q = self.real.carrier.wrap(q + v)
        raise ReductionError(f"after the fold at {self.drive} = {H:.6g} the state did not settle within "
                             f"t = {MAX_TIME:g}: it oscillates or escapes, which hysterons do not do")

    def step(self, S: np.ndarray, up: bool) -> Optional[Tuple[np.ndarray, float]]:
        k = int(S[0])
        key = (k, bool(up))
        if key in self._steps:
            return self._steps[key]
        b = self.branches[k]
        edge = b.hi if up else b.lo
        if not np.isfinite(edge):
            self._steps[key] = None
            return None
        H = edge + self.past if up else edge - self.past
        q_edge = b.q[-1] if up else b.q[0]
        u = self.rng.standard_normal(q_edge.shape)
        u *= self.perturb / np.linalg.norm(u)
        ends = [self._settle(q_edge, sgn * u, H) for sgn in (+1.0, -1.0)]
        if self.real.carrier.distance(ends[0][0], ends[1][0]) >= self.merge:
            raise ReductionError(
                f"state {k} loses stability at {self.drive} = {edge:.6g} in a way that the drive does not decide: "
                f"perturbations of size {self.perturb:.1g} in opposite directions lead to different states. The "
                "transition is a symmetric choice that noise makes (identical elements do this, through a pitchfork "
                "or a degenerate fold); a spread of their parameters lifts it")
        q_new, order = ends[0]
        t = self.identify(q_new, H)
        self.events.append({"state": k, "up": bool(up), "H": float(edge), "signs": b.signs.tolist(),
                            "order": order, "to": t})
        out = (np.array([float(t)]), float(H))
        self._steps[key] = out
        return out

    def interval(self, S: np.ndarray) -> Tuple[float, float]:
        b = self.branches[int(S[0])]
        return b.lo, b.hi


# ------------------------------------------------------------------------------------------------ the graph
def transition_graph(graph: EquationGraph) -> Dict[str, object]:
    """Every state reached from the low end of the range, with its up and down successors, in the form of
    counterpart.reachable (state 0 is the saturated start)."""
    for end in (graph.lo, graph.hi):
        found = graph.stable_states(end)
        if len(found) != 1:
            raise SpecError(f"slow_drive.range does not saturate the material: at {graph.drive} = {end:g} there are "
                            f"{len(found)} stable equilibria, not one; widen the range")
    if not graph.branches:
        graph._branch(graph.stable_states(graph.lo)[0], graph.lo)
    start = graph.identify(graph.stable_states(graph.lo)[0], graph.lo)
    order, index = [start], {start: 0}
    up: List[Optional[int]] = [None]
    down: List[Optional[int]] = [None]
    queue = deque([start])
    while queue:
        k = queue.popleft()
        for direction, edges in ((True, up), (False, down)):
            nxt = graph.step(np.array([float(k)]), direction)
            if nxt is None:
                continue
            t = int(nxt[0][0])
            if t not in index:
                index[t] = len(order)
                order.append(t)
                up.append(None)
                down.append(None)
                queue.append(t)
            edges[index[k]] = index[t]
    graph.m.n = len(graph.branches)
    # states are arrays [branch], so that counterpart.exact_return compares branches; the lists are indexed in the
    # order of discovery, with the saturated start first
    return {"states": [np.array([float(k)]) for k in order], "up": up, "down": down, "count": len(order),
            "graph": graph, "branches": order}


def survey(graph: EquationGraph, n_H: int = 24, n_starts: int = 32) -> Dict[str, object]:
    """Every stable equilibrium found from random starts on a grid of the drive, reachable or not, followed to its
    folds, and the fold events of all of them. The reduction needs the states that no history reaches: they carry
    the thresholds of the elements in configurations that the drive does not visit."""
    for H in np.linspace(graph.lo, graph.hi, n_H + 2)[1:-1]:
        for q in graph.stable_states(float(H), n_starts):
            graph.identify(q, float(H))
    unresolved = []
    k = 0
    while k < len(graph.branches):
        for up in (True, False):
            try:
                graph.step(np.array([float(k)]), up)
            except ReductionError as error:
                unresolved.append({"state": k, "up": up, "reason": str(error)})
        k += 1
    return {"states": len(graph.branches), "unresolved_folds": unresolved}


# ------------------------------------------------------------------------------------------------ reduction
def reduce(graph: EquationGraph) -> Dict[str, object]:
    """Fields, half-widths, couplings and drive signs of the hysteron model that reproduces the folds, by least
    squares, with the residual and the couplings that the folds do not determine."""
    n = graph.n_elements
    eta = np.zeros(n)
    rows: List[List[Tuple[np.ndarray, float]]] = [[] for _ in range(n)]
    hidden, mixed = 0, []
    for e in graph.events:
        if not e["order"]:
            hidden += 1
            continue
        i = e["order"][0]
        s = np.array(e["signs"], float)
        direction = -s[i] if e["up"] else s[i]           # the sign of eta_i that this event implies
        if eta[i] == 0:
            eta[i] = direction
        elif eta[i] != direction and i not in mixed:
            mixed.append(int(i))
        # unknowns of element i: h_i, b_i, then J_ij for j != i
        others = [j for j in range(n) if j != i]
        a = np.array([s[i], 1.0] + [s[i] * s[j] for j in others])
        rows[i].append((a, -s[i] * direction * e["H"]))
    h, b = np.zeros(n), np.zeros(n)
    J = np.zeros((n, n))
    determined = np.zeros((n, n), bool)
    residual, undetermined = [], []
    for i in range(n):
        others = [j for j in range(n) if j != i]
        if not rows[i]:
            undetermined.append({"element": i, "unknowns": "all: the element never switches first at a fold"})
            continue
        A = np.array([r[0] for r in rows[i]])
        y = np.array([r[1] for r in rows[i]])
        x, *_ = np.linalg.lstsq(A, y, rcond=None)
        residual += list(A @ x - y)
        h[i], b[i] = x[0], x[1]
        for c, j in enumerate(others):
            J[i, j] = x[2 + c]
        # an unknown is determined when its unit vector lies in the row space of A
        _, sv, vt = np.linalg.svd(A)
        rank = int(np.sum(sv > 1e-9 * max(1.0, sv.max())))
        basis = vt[:rank]
        for c in range(A.shape[1]):
            e_c = np.zeros(A.shape[1])
            e_c[c] = 1.0
            ok = np.linalg.norm(e_c - basis.T @ (basis @ e_c)) < 1e-6
            if c >= 2:
                determined[i, others[c - 2]] = ok
            elif not ok:
                undetermined.append({"element": i, "unknowns": "h and b" if c == 0 else "b"})
    for i in range(n):
        for j in range(n):
            if i != j and not determined[i, j]:
                undetermined.append({"element": i, "coupling_from": j})
    scale = float(np.mean(np.abs(b))) if np.any(b) else 1.0
    return {"fields": h.tolist(), "half_widths": b.tolist(), "couplings": J.tolist(), "drive": eta.tolist(),
            "determined": determined.tolist(), "undetermined": undetermined,
            "residual": float(np.sqrt(np.mean(np.square(residual)))) if residual else 0.0,
            "residual_relative": (float(np.sqrt(np.mean(np.square(residual)))) / scale) if residual else 0.0,
            "folds_that_switch_no_element": hidden, "elements_driven_both_ways": mixed}


def hysteron_model(reduced: Dict[str, object], name: str):
    """The reduced model as fieldbridge.memory.hysterons.Hysterons (reciprocal when J is symmetric)."""
    from .hysterons import Hysterons
    J = np.array(reduced["couplings"])
    n = J.shape[0]
    reciprocal = bool(np.allclose(J, J.T, atol=1e-3 * max(1.0, float(np.abs(J).max()))))
    if reciprocal:
        Js = 0.5 * (J + J.T)
        bonds = [(i, j, float(Js[i, j])) for i in range(n) for j in range(n) if i != j and Js[i, j] != 0]
    else:
        bonds = [(i, j, float(J[i, j])) for i in range(n) for j in range(n) if i != j and J[i, j] != 0]
    return Hysterons(name=f"{name} (reduced to hysterons)", n=n, bonds=bonds, h=np.array(reduced["fields"]),
                     b=np.maximum(np.array(reduced["half_widths"]), 0.0), eta=np.array(reduced["drive"]),
                     reciprocal=reciprocal, geometry={"shape": "graph"})


def check(spec: Dict, real: Realization, seed: int = 0, steps: int = STEPS, n_H: int = 24,
          n_starts: int = 32) -> Dict[str, object]:
    """The transition graph of the equations, its count of states and exact return; the reduced hysteron model,
    its structural prediction and its own count and return; and whether the two agree."""
    from . import counterpart as cp
    from .hysterons import predict
    drive, H_range, elements = from_spec(spec, real)
    graph = EquationGraph(real, drive, H_range, elements, steps=steps, seed=seed)
    tg = transition_graph(graph)
    depth = cp.reversals(tg)
    ret = cp.exact_return(tg)
    count = int(tg["count"])
    shallow = int(sum(d < cp.DEEP for d in depth))
    states = []
    for k, (branch, d) in enumerate(zip(tg["branches"], depth)):
        b = graph.branches[branch]
        states.append({"signs": "".join("+" if s > 0 else "-" for s in b.signs), "interval": [b.lo, b.hi],
                       "reversals": int(d), "up": tg["up"][k], "down": tg["down"][k], "continuous": b.continuous})
    sv = survey(graph, n_H=n_H, n_starts=n_starts)
    reached = set(tg["branches"])
    unreached = ["".join("+" if s > 0 else "-" for s in b.signs) for k, b in enumerate(graph.branches)
                 if k not in reached]
    out: Dict[str, object] = {
        "drive": drive, "range": list(H_range), "elements": elements, "states": states, "count": count,
        "C": float(np.log2(count)), "deep": float(np.log2(count) - np.log2(max(shallow, 1))),
        "max_reversals": int(max(depth)), "return": ret, "stable_states_not_reached": unreached,
        "unresolved_folds": sv["unresolved_folds"], "folds": len(graph.events)}
    continuous = sorted({i for b in graph.branches for i in b.continuous})
    if continuous:
        out["continuous_elements"] = continuous
    red = reduce(graph)
    out["reduced"] = red
    model = hysteron_model(red, str(spec.get("name", "material")))
    pred = predict(model)
    out["prediction"] = {k: pred[k] for k in ("guaranteed", "couplings_balanced", "statement")}
    out["frustration"] = {k: v for k, v in cp.frustration(model).items() if k in ("F", "exact")}
    measured = cp.measure(model)
    out["reduced_measure"] = {k: measured.get(k) for k in ("count", "return", "censored", "cycle")}
    out["agreement"] = {
        "count": measured.get("count") == count,
        "return": (measured.get("return") or {}).get("failed", -1) == 0 if ret["failed"] == 0
        else (measured.get("return") or {}).get("failed", 0) > 0}
    signs_known = not any("coupling_from" in u for u in red["undetermined"])
    if ret["failed"]:
        verdict = (f"{ret['failed']} of {ret['tested']} excursions do not return to their turning point in the "
                   "equations")
    else:
        verdict = f"every one of {ret['tested']} excursions returns to its turning point in the equations"
    if pred["guaranteed"] and signs_known:
        verdict += "; the reduced network is balanced with the drive, so the return is guaranteed"
        out["consistent"] = ret["failed"] == 0
    else:
        verdict += ("; the reduced network has a frustrated loop through the drive, which allows a failure"
                    if signs_known else "; some couplings are not determined by the folds, so the structural rule "
                    "is not decided")
        out["consistent"] = True
    out["verdict"] = verdict
    return out


def markdown(spec: Dict, r: Dict[str, object]) -> str:
    red = r["reduced"]
    lines = [f"# Return to a turning point of {spec.get('name', 'the material')}", "",
             str(spec.get("question", "")), "",
             f"The drive is {r['drive']} from {r['range'][0]:g} to {r['range'][1]:g}; the elements are "
             + ", ".join(f"sign({e})" for e in r["elements"]) + ".", "",
             f"- {r['verdict']}.",
             f"- States reached by some history of the drive: {r['count']} ({r['C']:.2f} bits); the most reversals "
             f"needed: {r['max_reversals']}.",
             f"- Reduced to hysterons: residual {red['residual_relative']:.3f} of the mean half-width; the reduced "
             f"model has {r['reduced_measure'].get('count')} states (equations: {r['count']}).",
             f"- Structure of the reduced network: {r['prediction']['statement']}", "",
             "| state | signs | interval of the drive | reversals | up to | down to |",
             "| --- | --- | --- | --- | --- | --- |"]
    for k, s in enumerate(r["states"]):
        lo, hi = s["interval"]
        span = (f"below {hi:.4g}" if not np.isfinite(lo) else f"above {lo:.4g}" if not np.isfinite(hi)
                else f"{lo:.4g} to {hi:.4g}")
        up = "none (end of the range)" if s["up"] is None else s["up"]
        down = "none (end of the range)" if s["down"] is None else s["down"]
        lines.append(f"| {k} | {s['signs']} | {span} | {s['reversals']} | {up} | {down} |")
    J = np.array(red["couplings"])
    lines += ["", "Reduced hysteron model (fields and couplings in units of the drive):", "",
              "| element | field h | half-width b | drive sign | couplings J from the others |",
              "| --- | --- | --- | --- | --- |"]
    for i, e in enumerate(r["elements"]):
        Jrow = ", ".join(f"{J[i, j]:+.4f}" + ("" if red["determined"][i][j] else " (not determined)")
                         for j in range(len(r["elements"])) if j != i)
        lines.append(f"| {e} | {red['fields'][i]:+.4f} | {red['half_widths'][i]:.4f} | {red['drive'][i]:+.0f} | "
                     f"{Jrow} |")
    lines += ["", "Statements about the noise-free, quasi-static dynamics of the supplied equations."]
    return "\n".join(lines) + "\n"
