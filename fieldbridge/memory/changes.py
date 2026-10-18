"""Changes of a model and their composition: routes from a published model to a mechanism it does not show.

Consecutive equations in papers differ by a few kinds of change: a term is added, removed or modified, a parameter
becomes the control, the range or the order of a protocol changes, a closure or an observable is replaced. Applied to
a supplied model, each change gives a new model whose mechanism the derivations of the memory language decide, so a
change is a derivation step that needs no target chosen in advance. Four kinds of change are implemented on equation
specifications:

  remove    drop one additive term of one equation;
  control   make another parameter the control, swept from a quarter to four times its value;
  range     extend the range of the control by ten times its width below it, above it, or on both sides;
  add       add to one equation one additive term of a donor model's equation in the same position, written in the
            model's own variables, with the donor's parameters set to one and a new coefficient g_add (only when a
            donor is given).

A changed model whose trajectories run away from large starts is discarded: the change removed the saturation and
there is no stable state to write. Every other changed model is judged by the symmetric-write and threshold-write
derivations (codiscovery), each in a child process of its own under a wall-clock limit: an exception raised inside a
solver's callback, or an alarm signal, would otherwise end the interpreter. The outcome of a model is the set of writes it reaches together with, for each write it
does not reach, the class of the reason why the derivation stops (no write point, a single stable state, an
oscillation, a one-sided write that needs a second parameter, a jump to a distant branch, several critical
directions, a vanishing cubic term). A change is a step when it reaches a write that the model did not reach or
changes one of these classes.

A pair of changes closes when its outcome is the outcome of the model or of one of its single changes, or combines
classes that single changes produce. It is open when it reaches a write that no single change reaches, or a class that
no single change produces: the two ordinary steps then do something that neither does alone. A derivation that fails
numerically or within its time limit decides nothing and does not open a pair.

A write reached after a change is a property of the changed model, not of the supplied one; every step names the
term or the parameter that was changed.
"""
from __future__ import annotations

import json
import re
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import codiscovery as cd
from . import spec as mspec

TARGETS = ("symmetric-write", "threshold-write")
KINDS = ("remove", "control", "range")
FUNCTIONS = ("exp", "sin", "cos", "tanh", "sqrt", "log", "abs", "cosh", "sinh")

# The class of the reason why a derivation stops, read from the obstruction text of codiscovery.
REASONS = (
    ("no write point along the control", "no write point"),
    ("no control parameter", "no write point"),
    ("a single stable state", "single stable state"),
    ("no stable state", "oscillation"),
    ("no oscillation", "no oscillation"),
    ("neutral cycles", "neutral cycles"),
    ("oscillat", "oscillation"),
    ("complex pair", "oscillation"),
    ("one-sided write", "one-sided write"),
    ("fold", "one-sided write"),
    ("positive cubic", "distant branch"),
    ("degenerate critical", "several critical directions"),
    ("directions become critical", "several critical directions"),
    ("vanishing", "vanishing cubic"),
    ("multiplies the whole drift", "the control scales the drift"),
    ("rate of time", "the control scales the drift"),
    ("time limit", "not derived"),
    ("LinAlgError", "not derived"),
    ("SVD", "not derived"),
)
UNDECIDED = ("not derived", "unclassified")


def reason_class(text: str) -> str:
    t = (text or "").lower()
    for key, kind in REASONS:
        if key.lower() in t:
            return kind
    return "unclassified"


# ------------------------------------------------------------------------------------------------ terms of a drift
def split_terms(expr: str) -> List[str]:
    """Top-level additive terms of an expression string, each with its sign."""
    terms, depth, cur = [], 0, ""
    for ch in expr:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if depth == 0 and ch in "+-" and cur.strip():
            prev = cur.rstrip()
            if prev[-1] not in "*/(^eE" or (prev[-1] in "eE" and not re.search(r"\d[eE]$", prev)):
                terms.append(cur.strip())
                cur = ch
                continue
        cur += ch
    if cur.strip():
        terms.append(cur.strip())
    return [t for t in terms if t not in ("+", "-")]


def join_terms(terms: Sequence[str]) -> str:
    out = ""
    for t in terms:
        t = t.strip()
        if not out:
            out = t[1:].strip() if t.startswith("+") else t
        elif t[0] in "+-":
            out += f" {t[0]} {t[1:].strip()}"
        else:
            out += " + " + t
    return out or "0"


# ------------------------------------------------------------------------------------------------ boundedness
def bounded(spec: Dict, radius: float = 4.0, directions: int = 32, horizon: float = 20.0, max_steps: int = 4000,
            seed: int = 0) -> Tuple[bool, str]:
    """Whether large states stay bounded, with the control at both ends of its range: short integrations from a ring
    of starts at `radius` times the carrier scale must neither run away to ten times that radius nor still be growing
    at the end. A change that removes the saturation runs away and leaves no stable state to write. The test
    integrates, because a drift averaged over directions misses a single run-away axis. Angles are bounded."""
    try:
        real = mspec.load(spec)
    except Exception as e:  # noqa: BLE001
        return False, f"not read: {type(e).__name__}"
    if getattr(real.carrier, "kind", "euclid") == "torus":
        return True, ""
    rng = np.random.default_rng(seed)
    dim = real.carrier.dim
    u = rng.normal(size=(directions, dim))
    u /= np.maximum(np.linalg.norm(u, axis=1, keepdims=True), 1e-12)
    if getattr(real.carrier, "kind", "") == "orthant":
        u = np.abs(u)
    scale = float(getattr(real.carrier, "scale", 1.0) or 1.0)
    ctrl = spec.get("control") if isinstance(spec.get("control"), dict) else None
    values: List[Optional[float]] = [None]
    if ctrl and ctrl.get("name") in real.params and isinstance(ctrl.get("range"), list) and len(ctrl["range"]) == 2:
        values = [float(ctrl["range"][0]), float(ctrl["range"][1])]
    R0 = radius * scale
    for v in values:
        where = f" for {ctrl['name']} = {v:.4g}" if v is not None else ""
        p = dict(real.params)
        if v is not None:
            p[ctrl["name"]] = v
        q = R0 * u.copy()
        t, steps, trail = 0.0, 0, []
        while t < horizon and steps < max_steps:
            with np.errstate(all="ignore"):
                f = real.drift(q, p)
            if not np.all(np.isfinite(f)):
                return False, f"unbounded: the drift is not finite at amplitude {R0:.3g}{where}"
            speed = float(np.max(np.linalg.norm(f, axis=1)))
            h = min(0.05, 0.05 * R0 / max(speed, 1e-12))
            with np.errstate(all="ignore"):
                k2 = real.drift(q + 0.5 * h * f, p)
                k3 = real.drift(q + 0.5 * h * k2, p)
                k4 = real.drift(q + h * k3, p)
            q = q + (h / 6.0) * (f + 2 * k2 + 2 * k3 + k4)
            t += h
            steps += 1
            norms = np.linalg.norm(q, axis=1)
            if not np.all(np.isfinite(norms)) or np.any(norms > 10 * R0):
                return False, f"unbounded: a state at amplitude {R0:.3g} runs away{where}"
            trail.append((t, float(np.max(norms))))
        if trail:
            end = trail[-1][1]
            earlier = next((n for tt, n in reversed(trail) if tt <= 0.75 * trail[-1][0]), trail[0][1])
            if end > 1.5 * R0 and end > 1.2 * earlier:
                return False, f"unbounded: states at amplitude {R0:.3g} are still growing at t = {trail[-1][0]:.3g}{where}"
    return True, ""


# ------------------------------------------------------------------------------------------------ single changes
def single_changes(spec: Dict, kinds: Sequence[str] = KINDS, rejected: Optional[List[Dict]] = None,
                   donor: Optional[Dict] = None, check_bounded: bool = True) -> List[Tuple[str, str, Dict]]:
    """(kind, description, changed specification) for every single change of an equation specification that keeps
    its states bounded; the discarded ones are appended to `rejected` with the reason. With a `donor`, the changes
    of kind add are included."""
    params = dict(spec.get("parameters") or {})
    drift = dict(spec.get("drift") or {})
    ctrl = spec.get("control") if isinstance(spec.get("control"), dict) else None
    out = []
    if "control" in kinds:
        for name, value in params.items():
            if ctrl and ctrl.get("name") == name:
                continue
            v = float(value)
            lo, hi = (v / 4, 4 * v) if v > 0 else ((4 * v, v / 4) if v < 0 else (-1.0, 1.0))
            s = dict(spec)
            s["control"] = {"name": name, "range": [lo, hi]}
            s["material"] = [p for p in params if p != name]
            out.append(("control", f"control = {name} over [{lo:.4g}, {hi:.4g}]", s))
    if "range" in kinds and ctrl and isinstance(ctrl.get("range"), list) and len(ctrl["range"]) == 2:
        lo, hi = float(ctrl["range"][0]), float(ctrl["range"][1])
        w = 10 * abs(hi - lo)
        for label, (a, b) in (("wider", (lo - w, hi + w)), ("above", (hi, hi + w)), ("below", (lo - w, lo))):
            s = dict(spec)
            s["control"] = {"name": ctrl["name"], "range": [a, b]}
            out.append(("range", f"range of {ctrl['name']} {label}: [{a:.4g}, {b:.4g}]", s))
    if "remove" in kinds:
        for var, expr in drift.items():
            terms = split_terms(str(expr))
            if len(terms) < 2:
                continue
            for i, t in enumerate(terms):
                s = dict(spec)
                s.pop("potential", None)   # a potential written for the old drift no longer belongs to it
                d = dict(drift)
                d[var] = join_terms(terms[:i] + terms[i + 1:])
                s["drift"] = d
                out.append(("remove", f"remove '{t}' from d{var}/dt", s))
    if donor is not None:
        out.extend(donor_terms(spec, donor))
    if not check_bounded:
        return out
    kept = []
    for kind, what, s in out:
        ok, why = bounded(s)
        if ok:
            kept.append((kind, what, s))
        elif rejected is not None:
            rejected.append({"kind": kind, "what": what, "why": why})
    return kept


def donor_terms(spec: Dict, donor: Dict) -> List[Tuple[str, str, Dict]]:
    """Changes of kind add: the first four additive terms of each equation of the donor, moved to the model's variable
    in the same position (the last one when the donor has more variables), the donor's parameters set to one and a
    new parameter g_add = 1 as the coefficient. A term that keeps a symbol which is neither a variable of the model nor
    a known function is skipped."""
    drift = dict(spec.get("drift") or {})
    params = dict(spec.get("parameters") or {})
    ctrl = spec.get("control") if isinstance(spec.get("control"), dict) else None
    fvars = list((donor.get("drift") or {}).keys())
    cvars = list(drift.keys())
    out: List[Tuple[str, str, Dict]] = []
    if not fvars or not cvars:
        return out
    mapping = {fv: cvars[min(i, len(cvars) - 1)] for i, fv in enumerate(fvars)}
    for fv in fvars:
        for t in split_terms(str(donor["drift"][fv]))[:4]:
            mapped = t.lstrip("+-").strip()
            for p in sorted((donor.get("parameters") or {}), key=len, reverse=True):
                mapped = re.sub(rf"(?<![A-Za-z0-9_]){re.escape(p)}(?![A-Za-z0-9_])", "1", mapped)
            for fvv in sorted(fvars, key=len, reverse=True):
                mapped = re.sub(rf"(?<![A-Za-z0-9_]){re.escape(fvv)}(?![A-Za-z0-9_])", mapping[fvv], mapped)
            words = re.findall(r"[A-Za-z_]\w*", mapped)
            if any(w not in cvars and w not in FUNCTIONS and w != "pi" for w in words):
                continue
            var = mapping[fv]
            s = dict(spec)
            s.pop("potential", None)
            d = dict(drift)
            d[var] = f"{drift[var]} + g_add*({mapped})"
            s["drift"] = d
            pr = dict(params)
            pr["g_add"] = 1.0
            s["parameters"] = pr
            s["material"] = [p for p in pr if not (ctrl and ctrl.get("name") == p)]
            out.append(("add", f"add g_add*({mapped}) to d{var}/dt (from {donor.get('name', 'the donor')})", s))
    return out


# ------------------------------------------------------------------------------------------------ judging a model
def _row(status: str, obstruction: str, seconds: float) -> Dict[str, object]:
    return {"status": status, "obstruction": obstruction[:200], "seconds": round(float(seconds), 1),
            "class": "reached" if status.startswith("reached") else reason_class(obstruction)}


def judge(spec: Dict, budget: int = 60, seed: int = 0, isolate: bool = True) -> Dict[str, Dict[str, object]]:
    """The symmetric-write and threshold-write derivations of a specification, without the law constants. With
    `isolate` (the default) each derivation runs in a child process of its own and is stopped after `budget` seconds
    of wall-clock time; without it the derivations run in this process with no time limit."""
    if not isolate:
        from . import codiscovery as cd
        import time
        try:
            real = mspec.load(spec)
        except Exception as e:  # noqa: BLE001
            return {t: _row("not read", f"{type(e).__name__}: {str(e)[:100]}", 0.0) for t in TARGETS}
        out = {}
        for t in TARGETS:
            t0 = time.time()
            try:
                row = cd.TARGETS[t](real, np.random.default_rng(seed), check_law=False)
                out[t] = _row(str(row.get("status", "")), str(row.get("obstruction", "")), time.time() - t0)
            except Exception as e:  # noqa: BLE001
                out[t] = _row("not derived", f"{type(e).__name__}: {str(e)[:120]}", time.time() - t0)
        return out
    from .cli import _child_step
    settings = {"seed": seed, "law": False, "trajectories": 0, "samples": 0, "trajectories_completion": 0}
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "model.json"
        path.write_text(json.dumps(spec), encoding="utf-8")
        for t in TARGETS:
            child, why = _child_step(path, settings, t, budget)
            if child is not None and child.get("not_read"):
                out[t] = _row("not read", str(child["not_read"]), 0.0)
            elif child is not None and t in child.get("targets", {}):
                r = child["targets"][t]
                out[t] = _row(str(r.get("status", "")), str(r.get("obstruction", r.get("note", ""))),
                              r.get("seconds", 0.0) or 0.0)
            else:
                out[t] = _row("not derived", why or "the derivation process returned nothing", float(budget))
    return out


def outcome(judged: Dict[str, Dict[str, object]]) -> Tuple[Tuple[str, ...], Tuple[Tuple[str, str], ...]]:
    """(writes reached, ((target, class), ...) for the others): the key by which outcomes are compared."""
    reached = tuple(sorted(t for t, r in judged.items() if r["class"] == "reached"))
    classes = tuple(sorted((t, str(r["class"])) for t, r in judged.items() if r["class"] != "reached"))
    return reached, classes


def closes(key, known: Sequence) -> Tuple[bool, str]:
    """Whether the outcome of a pair is contained in the outcomes of the model and of its single changes."""
    r_c, c_c = key
    if key in known:
        return True, "equals a known outcome"
    known_reached = {t for r, _ in known for t in r}
    new_targets = [t for t in r_c if t not in known_reached]
    if new_targets:
        return False, "reaches the " + " and the ".join(new_targets) + ", which no single change reaches"
    known_classes = {(t, k) for _, cl in known for t, k in cl}
    new_classes = [f"{t}: {k}" for t, k in c_c if (t, k) not in known_classes and k not in UNDECIDED]
    if new_classes:
        return False, "a reason for stopping that no single change produces (" + "; ".join(new_classes) + ")"
    if any(k in UNDECIDED for _, k in c_c):
        return True, "undecided: a derivation did not resolve"
    return True, "a combination of known outcomes"


# ------------------------------------------------------------------------------------------------ the routes
def changes(spec: Dict, budget: int = 60, seed: int = 0, pairs: bool = False, first: int = 8,
            kinds: Sequence[str] = KINDS, only: Optional[Sequence[str]] = None, progress=None,
            donor: Optional[Dict] = None, isolate: bool = True) -> Dict[str, object]:
    """Every single change of a specification, judged; with `pairs`, every pair formed from up to `first` first
    changes (steps first) followed by every change available on the changed model, two extensions of the same range
    excepted. `only` restricts the first changes to those whose description contains one of the given strings. With
    a `donor`, the single changes include its terms (kind add)."""
    base = judge(spec, budget, seed, isolate)
    key0 = outcome(base)
    rejected: List[Dict] = []
    singles = single_changes(spec, kinds, rejected, donor=donor)
    rows = []
    for kind, what, s in singles:
        j = judge(s, budget, seed, isolate)
        k = outcome(j)
        reached = sorted(set(k[0]) - set(key0[0]))
        changed = {t: [dict(key0[1]).get(t), dict(k[1]).get(t)] for t in dict(k[1])
                   if t in dict(key0[1]) and dict(k[1])[t] != dict(key0[1])[t]}
        rows.append({"kind": kind, "what": what, "judged": j, "reaches": reached, "changed": changed,
                     "step": bool(reached or changed), "spec": s})
        if progress:
            progress(f"{what}: {', '.join(reached) or ', '.join(f'{t} {a} -> {b}' for t, (a, b) in changed.items()) or 'no change'}")
    out: Dict[str, object] = {"name": spec.get("name"), "base": base, "singles": rows, "rejected": rejected,
                              "steps": sum(r["step"] for r in rows)}
    if not pairs:
        return out
    known = [key0] + [outcome(r["judged"]) for r in rows]
    firsts = [r for r in rows if r["step"]] + [r for r in rows if not r["step"]]
    if only:
        firsts = [r for r in rows if any(o in r["what"] for o in only)]
    firsts = firsts[:first]
    pair_rows = []
    for r1 in firsts:
        seconds = single_changes(r1["spec"], kinds)
        for kind, what, s in seconds:
            if r1["kind"] == "range" and kind == "range":
                continue
            j = judge(s, budget, seed, isolate)
            ok, why = closes(outcome(j), known)
            pair_rows.append({"first": [r1["kind"], r1["what"]], "second": [kind, what], "judged": j, "closes": ok,
                              "why": why, "spec": s})
            if progress and not ok:
                progress(f"open: {r1['what']} then {what}: {why}")
    out["pairs"] = pair_rows
    out["open"] = sum(not p["closes"] for p in pair_rows)
    return out


def markdown(name: str, res: Dict[str, object]) -> str:
    base = res["base"]
    lines = [f"# Changes of a model: {name}", "",
             "The supplied model: " + "; ".join(f"{t} {r['class']}" for t, r in base.items()) + ".", "",
             f"{len(res['singles'])} single changes keep the states bounded ({len(res['rejected'])} discarded because "
             f"states run away); {res['steps']} of them are steps.", "",
             "| change | symmetric write | threshold write | step |", "| --- | --- | --- | --- |"]
    for r in res["singles"]:
        j = r["judged"]
        lines.append(f"| {r['what']} | {j['symmetric-write']['class']} | {j['threshold-write']['class']} | "
                     f"{'reaches ' + ', '.join(r['reaches']) if r['reaches'] else ('yes' if r['step'] else '')} |")
    if "pairs" in res:
        P = res["pairs"]
        lines += ["", f"## Pairs of changes", "", f"{len(P)} pairs; {res['open']} open.", ""]
        for p in P:
            if not p["closes"]:
                lines.append(f"- {p['first'][1]}, then {p['second'][1]}: {p['why']}.")
    lines += ["", "A write reached after a change is a property of the changed model, not of the supplied one.", ""]
    return "\n".join(lines)
