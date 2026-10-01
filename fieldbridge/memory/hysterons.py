"""Return to a turning point: interacting hysterons under a slow drive (specification kind ``hysterons``).

Elements s_i = +1 or -1 switch at thresholds of their local field

    f_i = sum_j J_ij s_j + h_i + eta_i H,

and element i keeps its state while s_i f_i + b_i >= 0. Here h_i is a quenched random field, eta_i the sign and
strength with which the drive H acts on i, and b_i >= 0 half the width of the element's own loop: b = 0 is the
zero-temperature random-field Ising model (Sethna et al., Phys. Rev. Lett. 70, 3347, 1993), and J = 0 with b > 0 is
a set of independent hysterons (the Preisach model). The drive changes quasi-statically and without thermal noise.
Between two events nothing moves; at an event the element that has become unstable switches, and the avalanche that
follows is relaxed at fixed H, switching the least stable element first or, as a control, a randomly chosen unstable
one. The state is therefore a function of the sequence of turning points of the drive, not of its rate.

Return-point memory: when the drive turns at H1, makes an excursion that does not pass H1 and comes back, the
configuration at H1 is recovered exactly.

Prediction from structure. Count the drive as an element joined to every driven element i with the sign of eta_i. If
no loop of this network is frustrated, the relabeling sigma_i = g_i s_i (g_i = +1 or -1, g_i = sign(eta_i) on driven
elements) makes every coupling non-negative and the drive uniform. The relabeled dynamics has the no-passing property
(Middleton, Phys. Rev. Lett. 68, 670, 1992; Sethna et al., 1993): a configuration that lies above another at every
element stays above it under the same history of the drive. While the drive rises only upward switches occur, so an
avalanche ends with every element switched at most once and its final state does not depend on the order of the
switches. The return to a turning point is then exact. The same sign condition makes a system with an input monotone
(Angeli and Sontag, IEEE Trans. Automat. Control 48, 1684, 2003).

A frustrated loop through the drive allows the return to fail, also when no loop of couplings is frustrated: on a
bipartite lattice every loop of antiferromagnetic couplings is even, but in a uniform field every loop through the
drive is frustrated. It does not force a failure: random antiferromagnetic chains started from a large field return
exactly (Deutsch, Dhar and Narayan, Phys. Rev. Lett. 92, 227203, 2004). The structural condition is sufficient and
not necessary, and check() measures what the dynamics does.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .networks import sign_balance
from .spec import SCHEMA, SpecError

KIND = "hysterons"
SHAPES = ("square", "chain")
DRIVES = ("uniform", "staggered", "random")
MAX_UNITS = 20000


class AvalancheError(RuntimeError):
    """An avalanche did not end: possible only with couplings that are not reciprocal and not cooperative."""


@dataclass
class Hysterons:
    """Elements, couplings, random fields, half-widths of the elements' loops, and the action of the drive."""
    name: str
    n: int
    bonds: List[Tuple[int, int, float]]      # (i, j, J): j acts on i with J; reciprocal bonds appear in both directions
    h: np.ndarray
    b: np.ndarray
    eta: np.ndarray
    reciprocal: bool = True
    geometry: Dict[str, object] = field(default_factory=dict)
    out_idx: np.ndarray = field(init=False, repr=False)
    out_w: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        # padded lists of the elements on which each element acts; padding points to the element itself with weight 0
        lists: List[Dict[int, float]] = [dict() for _ in range(self.n)]
        for i, j, J in self.bonds:
            lists[j][i] = lists[j].get(i, 0.0) + float(J)
        k = max((len(d) for d in lists), default=0) or 1
        self.out_idx = np.tile(np.arange(self.n)[:, None], (1, k))
        self.out_w = np.zeros((self.n, k))
        for j, d in enumerate(lists):
            for c, (i, J) in enumerate(sorted(d.items())):
                self.out_idx[j, c], self.out_w[j, c] = i, J

    def local(self, s: np.ndarray, H: float) -> np.ndarray:
        f = self.h + self.eta * H
        np.add.at(f, self.out_idx.ravel(), (self.out_w * s[:, None]).ravel())
        return f

    def undirected(self) -> List[Tuple[int, int, float]]:
        """Each coupling once, as (i, j, J) with i < j for reciprocal couplings and as given otherwise."""
        if not self.reciprocal:
            return list(self.bonds)
        seen, out = set(), []
        for i, j, J in self.bonds:
            key = (min(i, j), max(i, j))
            if key not in seen:
                seen.add(key)
                out.append((key[0], key[1], J))
        return out

    def saturation(self) -> float:
        """A drive below which every driven element is in the state that the drive favours, whatever the others."""
        driven = self.eta != 0
        if not driven.any():
            raise SpecError("No element is driven: give the drive a nonzero sign on at least one element")
        reach = self._reach()
        return -float(np.max((reach[driven] + np.abs(self.h[driven]) + self.b[driven]) / np.abs(self.eta[driven]))) - 1.0

    def _reach(self) -> np.ndarray:
        r = np.zeros(self.n)
        np.add.at(r, self.out_idx.ravel(), np.abs(self.out_w).ravel())
        return r


# ------------------------------------------------------------------------------------------------ specification
def _values(item, n: int, rng: np.random.Generator, what: str) -> np.ndarray:
    if item is None:
        return np.zeros(n)
    if isinstance(item, (int, float)):
        return np.full(n, float(item))
    if isinstance(item, list):
        if len(item) != n:
            raise SpecError(f"{what}: expected {n} values, one per element, got {len(item)}")
        return np.array([float(v) for v in item])
    if isinstance(item, dict):
        dist = item.get("distribution")
        if dist == "gaussian":
            return float(item.get("mean", 0.0)) + float(item["width"]) * rng.standard_normal(n)
        if dist == "uniform":
            lo, hi = float(item["low"]), float(item["high"])
            if hi < lo:
                raise SpecError(f"{what}: 'high' must not be below 'low'")
            return rng.uniform(lo, hi, n)
        raise SpecError(f"{what}: distribution must be 'gaussian' (width, mean) or 'uniform' (low, high)")
    raise SpecError(f"{what}: give a number, a list with one value per element, or a distribution")


def _lattice(lat: Dict, rng: np.random.Generator,
             couplings: Dict) -> Tuple[int, List[Tuple[int, int, float]], List[int], Dict]:
    shape, L = lat.get("shape"), int(lat.get("L", 0))
    if shape not in SHAPES:
        raise SpecError(f"lattice.shape must be one of {SHAPES}")
    if L < 3:
        raise SpecError("lattice.L must be at least 3")
    if shape == "square":
        n = L * L
        idx = np.arange(n).reshape(L, L)
        pairs = [(int(a), int(c)) for a, c in zip(idx.ravel(), np.roll(idx, -1, 1).ravel())]
        pairs += [(int(a), int(c)) for a, c in zip(idx.ravel(), np.roll(idx, -1, 0).ravel())]
        parity = [int((i // L + i % L) % 2) for i in range(n)]
    else:
        n = L
        pairs = [(i, (i + 1) % L) for i in range(L)]
        parity = [i % 2 for i in range(n)]
    value = float(couplings.get("value", 1.0))
    p = float(couplings.get("negative_fraction", 0.0))
    if not 0.0 <= p <= 1.0:
        raise SpecError("couplings.negative_fraction must lie between 0 and 1")
    signs = np.where(rng.random(len(pairs)) < p, -1.0, 1.0)
    bonds = [(i, j, value * s) for (i, j), s in zip(pairs, signs)]
    geometry = {"shape": shape, "L": L, "bipartite": L % 2 == 0}
    return n, bonds, parity, geometry


def from_spec(spec: Dict) -> Hysterons:
    """The model of a specification of kind ``hysterons``; SpecError says what to change."""
    if spec.get("schema") != SCHEMA:
        raise SpecError(f"schema must be {SCHEMA!r}")
    if spec.get("kind") != KIND:
        raise SpecError(f"kind must be {KIND!r} for this command")
    hs = spec.get("hysterons")
    if not isinstance(hs, dict):
        raise SpecError("A specification of kind 'hysterons' needs a 'hysterons' object")
    rng = np.random.default_rng(int(hs.get("seed", 0)))
    parity = None
    reciprocal = True
    if "lattice" in hs:
        n, pairs, parity, geometry = _lattice(hs["lattice"], rng, hs.get("couplings") or {})
        bonds = [(i, j, J) for i, j, J in pairs] + [(j, i, J) for i, j, J in pairs]
    elif "units" in hs:
        n = int(hs["units"])
        if n < 1:
            raise SpecError("units must be a positive number of elements")
        geometry = {"shape": "graph"}
        bonds = []
        reciprocal = not bool(hs.get("directed", False))
        for c in hs.get("couplings") or []:
            if isinstance(c, dict):
                i, j, J = int(c["to"]), int(c["from"]), float(c["J"])
            else:
                i, j, J = int(c[0]), int(c[1]), float(c[2])
            if not (0 <= i < n and 0 <= j < n):
                raise SpecError(f"coupling {c}: elements are numbered 0 to {n - 1}")
            if i == j:
                raise SpecError(f"coupling {c}: an element does not act on itself; use half_widths for its own loop")
            bonds += [(i, j, J)] + ([(j, i, J)] if reciprocal else [])
    else:
        raise SpecError("Give either 'lattice' (shape, L) or 'units' with a list of 'couplings'")
    if n > MAX_UNITS:
        raise SpecError(f"At most {MAX_UNITS} elements are supported")
    h = _values(hs.get("fields"), n, rng, "fields")
    b = _values(hs.get("half_widths"), n, rng, "half_widths")
    if (b < 0).any():
        raise SpecError("half_widths must not be negative: b is half the width of an element's own loop")
    drive = hs.get("drive", "uniform")
    if isinstance(drive, list):
        eta = _values(drive, n, rng, "drive")
    elif drive == "uniform":
        eta = np.ones(n)
    elif drive == "staggered":
        if parity is None:
            raise SpecError("A staggered drive needs a lattice; for a graph, give the drive as a list")
        eta = np.where(np.array(parity) == 0, 1.0, -1.0)
    elif drive == "random":
        eta = np.where(rng.random(n) < 0.5, 1.0, -1.0)
    else:
        raise SpecError(f"drive must be one of {DRIVES} or a list with one value per element")
    return Hysterons(name=str(spec.get("name", spec.get("id", "hysterons"))), n=n, bonds=bonds, h=h, b=b, eta=eta,
                     reciprocal=reciprocal, geometry=geometry)


# ------------------------------------------------------------------------------------------------ dynamics
class Run:
    """One history of the drive, started at saturation (a large negative drive)."""

    def __init__(self, model: Hysterons, order: str = "least", rng: Optional[np.random.Generator] = None):
        if order not in ("least", "random"):
            raise ValueError("order must be 'least' or 'random'")
        self.m, self.order = model, order
        self.rng = rng if rng is not None else np.random.default_rng(0)
        self.H = model.saturation()
        self.s = np.where(model.eta > 0, -1.0, np.where(model.eta < 0, 1.0, -1.0))
        self.f = model.local(self.s, self.H)
        self.switches = 0
        self.trace: Optional[List[Tuple[float, float]]] = None    # (H, response) after each event, when enabled
        self.events: Optional[List[Tuple[float, int]]] = None     # (H, switches so far) after each event, when enabled
        self.reference: Optional[np.ndarray] = None               # a configuration to compare with, when enabled
        self.unlike: List[int] = []                               # elements unlike the reference after each event
        self._relax()

    def clone(self, rng: Optional[np.random.Generator] = None) -> "Run":
        """The same configuration and drive, continued independently (with its own random order, if given)."""
        other = Run.__new__(Run)
        other.m, other.order, other.rng = self.m, self.order, rng if rng is not None else self.rng
        other.H, other.s, other.f, other.switches = self.H, self.s.copy(), self.f.copy(), self.switches
        other.trace = other.events = other.reference = None
        other.unlike = []
        return other

    def response(self) -> float:
        """The observable conjugate to the drive: sum_i eta_i s_i / sum_i |eta_i|."""
        return float(np.dot(self.m.eta, self.s) / np.abs(self.m.eta).sum())

    def _flip(self, i: int) -> None:
        self.s[i] = -self.s[i]
        self.f[self.m.out_idx[i]] += 2.0 * self.m.out_w[i] * self.s[i]
        self.switches += 1

    def _relax(self) -> None:
        cap, count = 50 * self.m.n + 1000, 0
        while True:
            st = self.s * self.f + self.m.b
            if self.order == "least":
                i = int(np.argmin(st))
                if st[i] >= 0:
                    return
            else:
                unstable = np.flatnonzero(st < 0)
                if unstable.size == 0:
                    return
                i = int(unstable[self.rng.integers(unstable.size)])
            self._flip(i)
            count += 1
            if count > cap:
                raise AvalancheError(f"an avalanche did not end after {cap} switches; with couplings that are not "
                                     "reciprocal the dynamics can cycle")

    def ramp(self, H_to: float) -> None:
        """Move the drive to H_to, switching at each event the element that becomes unstable and relaxing the
        avalanche that follows."""
        eta = self.m.eta
        direction = 1.0 if H_to > self.H else -1.0
        absn = np.abs(eta)
        while True:
            # elements whose stability the drive lowers: s_i eta_i direction < 0; they lose it after |dH| = st / |eta|
            cand = self.s * eta * direction < 0
            nxt = H_to
            if cand.any():
                st = self.s[cand] * self.f[cand] + self.m.b[cand]
                step = float(np.min(st / absn[cand]))
                if step < abs(H_to - self.H):
                    nxt = self.H + direction * max(step, 0.0)
            self.f += eta * (nxt - self.H)
            self.H = nxt
            if nxt == H_to:
                self._relax()
                self._record()
                return
            # at the event the threshold element is marginal; move past it by a negligible field
            self.f += eta * direction * 1e-12
            self.H += direction * 1e-12
            self._relax()
            self._record()

    def _record(self) -> None:
        if self.trace is not None:
            self.trace.append((self.H, self.response()))
        if self.events is not None:
            self.events.append((self.H, self.switches))
        if self.reference is not None:
            self.unlike.append(int(np.sum(self.s != self.reference)))


def major_loop(model: Hysterons) -> Dict[str, List[Tuple[float, float]]]:
    """The rising and falling branches between saturations, as (H, response) after each event."""
    run = Run(model)
    top = -run.H
    run.trace = [(run.H, run.response())]
    run.ramp(top)
    up = run.trace
    run.trace = [(run.H, run.response())]
    run.ramp(-top)
    return {"up": up, "down": run.trace}


def switching_window(model: Hysterons, lo: float = 0.1, hi: float = 0.9, falling: bool = False) -> Tuple[float, float]:
    """The drives between which a branch of the major loop makes the fractions lo to hi of its switches: the part of
    the loop in which avalanches occur. The rising branch starts at saturation; the falling branch (falling=True)
    starts at the top. Returned as (lower, upper)."""
    run = Run(model)
    bottom = run.H
    if falling:
        run.ramp(-bottom)
    start = run.switches
    run.events = []
    run.ramp(bottom if falling else -bottom)
    events = [(H, c - start) for H, c in run.events if c > start]
    if not events:
        raise SpecError("No element switches between the two saturations: the drive does not reach the thresholds")
    total = events[-1][1]
    fields = [H for H, _ in events]
    counts = np.array([c for _, c in events]) / total
    first = fields[int(np.searchsorted(counts, lo))]
    last = fields[min(int(np.searchsorted(counts, hi)), len(fields) - 1)]
    return float(min(first, last)), float(max(first, last))


def excursion_depth(H1: float, rising: Tuple[float, float], falling: Tuple[float, float]) -> Tuple[float, float]:
    """The range of the lower turning point H2 of an excursion from H1: from the field at which the falling branch has
    made 90 per cent of its switches up to the field of its first 10 per cent, and in any case below H1, so that the
    descent reaches the part of the loop in which elements switch back."""
    step = 0.05 * max(rising[1] - rising[0], 1e-6)
    upper = min(falling[1], H1 - step)
    return min(falling[0], upper - step), upper


# ------------------------------------------------------------------------------------------------ structure
def predict(model: Hysterons) -> Dict[str, object]:
    """Whether the return to a turning point is guaranteed, from the signs of the couplings and of the drive."""
    n = model.n
    couplings = [(i, j, int(np.sign(J))) for i, j, J in model.undirected() if J != 0]
    driven = np.flatnonzero(model.eta != 0)
    with_drive = couplings + [(int(i), n, int(np.sign(model.eta[i]))) for i in driven]
    extended = sign_balance(n + 1, with_drive)
    alone = sign_balance(n, couplings)
    g_drive = extended["gauge"][n]
    relabeling = [int(c * g_drive) for c in extended["gauge"][:n]]
    # the shortest loops through the drive: the drive, i and j for every coupling between driven elements
    pairs = [(i, j, s) for i, j, s in couplings if model.eta[i] != 0 and model.eta[j] != 0]
    frustrated = [(i, j) for i, j, s in pairs if s * np.sign(model.eta[i]) * np.sign(model.eta[j]) < 0]
    out: Dict[str, object] = {
        "guaranteed": bool(extended["balanced"]), "couplings_balanced": bool(alone["balanced"]),
        "reciprocal": model.reciprocal, "bonds": len(pairs), "bonds_frustrated_with_drive": len(frustrated),
        "frustrated_with_drive": len(frustrated) / len(pairs) if pairs else 0.0,
        "frustrated_plaquettes": frustrated_plaquettes(model)}
    if extended["balanced"]:
        out["relabeling"] = relabeling
        out["statement"] = ("No loop through the drive is frustrated. Relabeling the states as σᵢ = gᵢsᵢ makes every "
                            "coupling non-negative and the drive uniform; the dynamics then has no passing, and every "
                            "return to a turning point of the drive is exact.")
        out["canonical"] = canonical_form(model, relabeling)
    else:
        where = ("although no loop of couplings is frustrated" if alone["balanced"]
                 else "and loops of couplings are frustrated as well")
        out["statement"] = (f"{len(frustrated)} of {len(pairs)} couplings close a frustrated loop through the drive, "
                            f"{where}. The return to a turning point is not guaranteed: it can fail, and it need not.")
    return out


def frustrated_plaquettes(model: Hysterons) -> Optional[float]:
    """Fraction of the elementary squares of a square lattice whose product of coupling signs is negative."""
    if model.geometry.get("shape") != "square":
        return None
    L = int(model.geometry["L"])
    J = {(i, j): Jv for i, j, Jv in model.bonds}
    neg = 0
    for i in range(L * L):
        r, c = divmod(i, L)
        a, b_ = i, r * L + (c + 1) % L
        d = ((r + 1) % L) * L + c
        e = ((r + 1) % L) * L + (c + 1) % L
        prod = J[(a, b_)] * J[(b_, e)] * J[(e, d)] * J[(d, a)]
        neg += prod < 0
    return neg / (L * L)


def canonical_form(model: Hysterons, relabeling: Sequence[int]) -> Dict[str, object]:
    """The relabeled model: cooperative couplings g_i J_ij g_j >= 0 and the drive |eta_i| H, with fields g_i h_i."""
    g = np.asarray(relabeling, dtype=float)
    K = [g[i] * J * g[j] for i, j, J in model.undirected()]
    kind = ("independent hysterons (Preisach)" if not K or max(abs(k) for k in K) == 0
            else "cooperative hysterons in a uniform drive" + (" (random-field Ising ferromagnet)"
                                                              if np.all(model.b == 0) else ""))
    return {"form": kind, "min_coupling": float(min(K)) if K else 0.0,
            "drive_positive": bool(np.all(g * model.eta >= 0))}


def relabel(model: Hysterons, relabeling: Sequence[int]) -> Hysterons:
    """The same material in the variables sigma_i = g_i s_i."""
    g = np.asarray(relabeling, dtype=float)
    bonds = [(i, j, g[i] * J * g[j]) for i, j, J in model.bonds]
    return replace(model, bonds=bonds, h=g * model.h, eta=g * model.eta)


# ------------------------------------------------------------------------------------------------ measurement
def subloops(model: Hysterons, rng: np.random.Generator, count: int = 24, nested: bool = False,
             order: str = "least", window: Optional[Tuple[float, float]] = None,
             falling: Optional[Tuple[float, float]] = None) -> List[Dict[str, object]]:
    """Excursions H1 -> H2 -> H1 from the rising branch, with H1 inside the switching window of the rising branch and
    H2 inside that of the falling branch (excursion_depth); nested=True adds an inner excursion H2 -> H3 -> H2 with
    H2 < H3 < H1."""
    lo, hi = window if window is not None else switching_window(model)
    falling = falling if falling is not None else switching_window(model, falling=True)
    plan = []
    for _ in range(count):
        H1 = float(rng.uniform(lo, hi))
        H2 = float(rng.uniform(*excursion_depth(H1, (lo, hi), falling)))
        plan.append((H1, H2, H2 + float(rng.uniform(0.1, 0.9)) * (H1 - H2) if nested else None))
    seeds = rng.integers(1 << 31, size=count + 1)
    # every excursion starts from the rising branch: one rise from saturation, stopped at each H1 in turn
    rise = Run(model, order, np.random.default_rng(int(seeds[-1])))
    starts: Dict[int, Optional[Run]] = {}
    for k in sorted(range(count), key=lambda k: plan[k][0]):
        try:
            rise.ramp(plan[k][0])
            starts[k] = rise.clone(np.random.default_rng(int(seeds[k])))
        except AvalancheError:
            starts.update({j: None for j in range(count) if j not in starts})
            break
    out = []
    for k, (H1, H2, H3) in enumerate(plan):
        run = starts[k]
        if run is None:
            out.append({"H1": H1, "H2": H2, "H3": H3, "differ": None, "ended": False})
            continue
        s1 = run.s.copy()
        try:
            run.ramp(H2)
            if nested:
                run.ramp(H3)
                run.ramp(H2)
            run.ramp(H1)
        except AvalancheError:
            out.append({"H1": H1, "H2": H2, "H3": H3, "differ": None, "ended": False})
            continue
        out.append({"H1": H1, "H2": H2, "H3": H3, "differ": int(np.sum(run.s != s1)), "ended": True})
    return out


def excursion(model: Hysterons, H1: float, H2: float) -> Dict[str, object]:
    """One excursion H1 -> H2 -> H1 from the rising branch: the response after each event, the progress along the
    excursion (0 at H1, 1 at H2, 2 back at H1) and the number of elements unlike the state at H1."""
    run = Run(model)
    run.ramp(H1)
    s1 = run.s.copy()
    run.trace, run.reference, run.unlike = [(run.H, run.response())], s1, [0]
    run.ramp(H2)
    down = len(run.trace)
    run.ramp(H1)
    span = H1 - H2
    progress = ([(H1 - H) / span for H, _ in run.trace[:down]] + [1 + (H - H2) / span for H, _ in run.trace[down:]])
    return {"H1": H1, "H2": H2, "trace": run.trace, "progress": progress, "unlike": run.unlike,
            "differ": int(np.sum(run.s != s1))}


def check(model: Hysterons, rng: np.random.Generator, count: int = 24) -> Dict[str, object]:
    """The prediction and the measured returns: simple and nested subloops, least stable element first and in a
    random order."""
    pred = predict(model)
    try:
        window, falling = switching_window(model), switching_window(model, falling=True)
    except AvalancheError:
        return {"prediction": pred, "subloops": {}, "total": 0, "failed": 0, "did_not_end": 1,
                "verdict": "an avalanche on the rising branch did not end: without a rule for cycles the history of "
                           "the drive does not fix the state",
                "consistent": not pred["guaranteed"]}
    rows = {}
    for nested in (False, True):
        for order in ("least", "random"):
            key = f"{order}{'_nested' if nested else ''}"
            loops = subloops(model, rng, count, nested=nested, order=order, window=window, falling=falling)
            ended = [r for r in loops if r["ended"]]
            rows[key] = {"subloops": len(loops), "failed": sum(r["differ"] > 0 for r in ended),
                         "did_not_end": len(loops) - len(ended),
                         "max_differ": max((r["differ"] for r in ended), default=0)}
    total = sum(r["subloops"] for r in rows.values())
    failed = sum(r["failed"] for r in rows.values())
    endless = sum(r["did_not_end"] for r in rows.values())
    endless_text = f"; in {endless} an avalanche did not end" if endless else ""
    if pred["guaranteed"]:
        ok = failed == 0 and endless == 0
        verdict = ("consistent: every subloop returned, as predicted" if ok
                   else f"INCONSISTENT: {failed} of {total} subloops did not return{endless_text}, although the return "
                        "is guaranteed")
    else:
        ok = True
        verdict = (f"{failed} of {total} subloops did not return{endless_text}; the prediction allows this"
                   if failed or endless else
                   f"every one of {total} subloops returned: a frustrated loop through the drive allows a failure "
                   "without forcing one")
    return {"prediction": pred, "subloops": rows, "total": total, "failed": failed, "did_not_end": endless,
            "verdict": verdict, "consistent": ok}


# ------------------------------------------------------------------------------------------------ the card
def shown_excursion(model: Hysterons, rising: Tuple[float, float], falling: Tuple[float, float],
                    failing: bool = False) -> Dict[str, object]:
    """The excursion a figure shows: from the middle of the rising window into the middle of the falling one, or,
    with failing=True, the first of a grid of excursions that does not return (the middle one if none fails)."""
    lo, hi = rising

    def at(a: float, depth: float) -> Dict[str, object]:
        H1 = lo + a * (hi - lo)
        low, up = excursion_depth(H1, rising, falling)
        return excursion(model, H1, up - depth * (up - low))

    shown = at(0.6, 0.5)
    if failing:
        for a in (0.3, 0.45, 0.6, 0.75, 0.9):
            for depth in (0.3, 0.5, 0.7, 0.9):
                trial = at(a, depth)
                if trial["differ"]:
                    return trial
    return shown


def structure(model: Hysterons, patch: int = 8) -> Dict[str, object]:
    """The signs of the network with the drive, for a figure: a patch of the lattice with the sign of each element's
    drive and the sign of eta_i J_ij eta_j on each bond, or, without couplings, the switching fields of each element
    (its Preisach pair)."""
    shape = model.geometry.get("shape")
    if shape in ("square", "chain"):
        L = int(model.geometry["L"])
        side = min(L, patch) if shape == "square" else min(L, 3 * patch)
        pos = ({r * L + c: (c, r) for r in range(side) for c in range(side)} if shape == "square"
               else {i: (i, 0) for i in range(side)})
        bonds = [[*pos[i], *pos[j], int(np.sign(model.eta[i] * J * model.eta[j]))]
                 for i, j, J in model.undirected() if i in pos and j in pos and J != 0
                 and abs(pos[i][0] - pos[j][0]) + abs(pos[i][1] - pos[j][1]) == 1]
        return {"shape": shape, "nodes": [[x, y, int(np.sign(model.eta[i]))] for i, (x, y) in pos.items()],
                "bonds": bonds}
    if not model.bonds:
        driven = model.eta != 0
        up = (model.b - model.h) / np.where(driven, np.abs(model.eta), 1)
        down = -(model.b + model.h) / np.where(driven, np.abs(model.eta), 1)
        return {"shape": "independent", "preisach": [[float(d), float(u)] for d, u in zip(down[driven], up[driven])]}
    return {"shape": "graph"}


def card(model: Hysterons, rng: np.random.Generator, count: int = 24) -> Dict[str, object]:
    """Everything a figure or the page shows of one material: the prediction and the measured subloops (check, with
    this generator first), the class, the major loop, one excursion and the signs of the network with the drive."""
    result = check(model, rng, count)
    cls = "return-point" if result["failed"] == 0 and result["did_not_end"] == 0 else "no-return"
    out = {"kind": KIND, "name": model.name, "class": cls, "elements": model.n, "check": result,
           "structure": structure(model), "window": None, "falling": None, "loop": None, "excursion": None}
    try:
        rising, falling = switching_window(model), switching_window(model, falling=True)
        out.update(window=list(rising), falling=list(falling), loop=major_loop(model),
                   excursion=shown_excursion(model, rising, falling, failing=cls == "no-return"))
    except AvalancheError:
        pass
    return out
