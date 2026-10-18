"""Typed moves on a drift specification, and the closure of their compositions.

A move changes one slot of a specification of kind `equations`: the protocol (which parameter is the control and over
which range) or the operator apparatus (one additive term of one equation dropped, or one term added from a donor
specification). The four moves, with the transition types of the literature they correspond to:

  add_protocol_order        declare a parameter as the control, swept over [value/4, 4 value];
  modify_protocol_order     widen or shift the control's range by a decade on either side;
  project_operator_apparatus    drop one additive term of one equation;
  add_operator_apparatus    add to one equation one additive term of the donor's equation in the same position, in the
                            specification's own variables, with a new parameter g_add for its coefficient.

A candidate is kept only when the modified realization is bounded: short integrations from a ring of starts at four
times the carrier scale, with the control at both ends of its range, must neither run away to ten times that amplitude
nor still be growing at the end. A change that removes the saturation (a loss rate made negative, the saturating term
dropped) has no stable state to write and is rejected; a static radial test fails both ways (a Kerr oscillator above
threshold has outward radial drift and stays bounded by its rotation; a frustrated triangle without one cubic has an
inward direction-averaged drift and runs away along one axis), so the test integrates.

Each kept candidate is judged by the write targets of `codiscovery`: the symmetric and the threshold write, each
reached or obstructed, and the class of the obstruction (what a completion would have to keep). A step is a move that
reaches a target the specification did not reach, or changes an obstruction class. A composition of two moves closes
when its outcome equals the outcome of the specification or of a single step, or combines obstruction classes the
single steps already show; it is open when it reaches a target no single step reaches or lands in an obstruction class
no single step shows. Open compositions get a third move. Nothing here establishes that a modified model describes a
physical system; the report says what the rules find on it.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from . import spec as mspec

MOVES = ("add_protocol_order", "modify_protocol_order", "project_operator_apparatus", "add_operator_apparatus")
TARGETS = ("symmetric-write", "threshold-write")
# obstruction texts of the derivations, in the order tested, with what a completion keeps
OBSTRUCTION_CLASSES = [
    ("no write point along the control", ("none", "no write point: the control never destabilizes the state")),
    ("no control parameter", ("none", "no control parameter")),
    ("a single stable state", ("none", "single stable state: nothing to choose")),
    ("no stable state", ("phase", "the state oscillates: the phase-locking target applies")),
    ("oscillat", ("phase", "the state oscillates: the phase-locking target applies")),
    ("complex pair", ("phase", "a complex pair crosses: keep the phase")),
    ("no oscillation", ("none", "no oscillation: phase locking does not apply")),
    ("neutral cycles", ("none", "a family of neutral cycles")),
    ("one-sided write", ("second parameter", "one-sided write: keep a second parameter, sweep through the cusp")),
    ("fold", ("second parameter", "fold: keep a second parameter, sweep through the cusp")),
    ("positive cubic", ("distant branch", "positive cubic: keep the distant branch")),
    ("degenerate critical", ("critical space", "degenerate critical mode: keep the critical space")),
    ("directions become critical", ("critical space", "degenerate critical mode: keep the critical space")),
    ("vanishing", ("leading order", "vanishing cubic: keep the leading term of higher order")),
    ("multiplies the whole drift", ("none", "the control scales time: no write")),
    ("time limit", ("not derived", "time limit")),
    ("not derived", ("not derived", "not derived")),
    ("LinAlgError", ("not derived", "numerical failure of the reduction")),
    ("SVD", ("not derived", "numerical failure of the reduction")),
]
LAW_COMPLETIONS = {"second parameter", "distant branch", "critical space", "leading order"}
UNDECIDED = ("not derived", "unclassified")
FUNCTIONS = ("exp", "sin", "cos", "tanh", "sqrt", "log", "abs", "cosh", "sinh")


def obstruction_class(text: str) -> Dict[str, object]:
    """The class of an obstruction from the derivation's text: what a completion keeps, and whether it has a law."""
    t = (text or "").lower()
    for key, (kind, keeps) in OBSTRUCTION_CLASSES:
        if key.lower() in t:
            return {"kind": kind, "keeps": keeps, "has_law": kind in LAW_COMPLETIONS}
    return {"kind": "unclassified", "keeps": (text or "")[:120], "has_law": False}


# ------------------------------------------------------------------------------------------------ terms
def split_terms(expr: str) -> List[str]:
    """Top-level additive terms of an expression, each with its sign."""
    terms, depth, cur = [], 0, ""
    for ch in expr:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if depth == 0 and ch in "+-" and cur.strip():
            prev = cur.rstrip()
            if prev and prev[-1] not in "*/(^" and not re.search(r"[eE]$", prev):
                terms.append(cur.strip())
                cur = ch
                continue
        cur += ch
    if cur.strip():
        terms.append(cur.strip())
    return [t for t in terms if t not in ("+", "-")]


def join_terms(terms: List[str]) -> str:
    out = ""
    for t in terms:
        t = t.strip()
        if not out:
            out = t
        elif t[0] in "+-":
            out += " " + t[0] + " " + t[1:].strip()
        else:
            out += " + " + t
    return out or "0"


# ------------------------------------------------------------------------------------------------ boundedness
def bounded(spec: Dict, radius: float = 4.0, directions: int = 32, horizon: float = 20.0, max_steps: int = 4000,
            seed: int = 0) -> Tuple[bool, str]:
    """Whether the drift keeps large-amplitude states bounded, with the control at both ends of its range (see the
    module docstring). Angles (torus carriers) are bounded by construction."""
    try:
        real = mspec.load(spec)
    except Exception as e:  # noqa: BLE001
        return False, f"not read: {type(e).__name__}: {str(e)[:80]}"
    if real.carrier.kind == "torus":
        return True, ""
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(directions, real.carrier.dim))
    u /= np.maximum(np.linalg.norm(u, axis=1, keepdims=True), 1e-12)
    if real.carrier.kind == "orthant":
        u = np.abs(u)
    scale = float(real.carrier.scale or 1.0)
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
                return False, f"unbounded: a state at amplitude {R0:.3g} runs away{where} (the change removed the saturation)"
            trail.append((t, float(np.max(norms))))
        if trail:
            end = trail[-1][1]
            earlier = next((n for tt, n in reversed(trail) if tt <= 0.75 * trail[-1][0]), trail[0][1])
            if end > 1.5 * R0 and end > 1.2 * earlier:
                return False, f"unbounded: states at amplitude {R0:.3g} are still growing at t = {trail[-1][0]:.3g}{where}"
    return True, ""


# ------------------------------------------------------------------------------------------------ moves
def candidates(spec: Dict, donor: Optional[Dict] = None, rejected: Optional[List[Dict]] = None,
               check_bounded: bool = True) -> List[Dict]:
    """Every single move on `spec`, as {"move", "what", "spec"}; the candidates that the boundedness test rejects go to
    `rejected` with the reason."""
    params = dict(spec.get("parameters") or {})
    drift = dict(spec.get("drift") or {})
    ctrl = spec.get("control") if isinstance(spec.get("control"), dict) else None
    out: List[Dict] = []
    for name, value in params.items():
        if ctrl and ctrl.get("name") == name:
            continue
        v = float(value)
        lo, hi = (v / 4, 4 * v) if v > 0 else ((4 * v, v / 4) if v < 0 else (-1.0, 1.0))
        s = dict(spec)
        s["control"] = {"name": name, "range": [lo, hi]}
        s["material"] = [p for p in params if p != name]
        out.append({"move": "add_protocol_order", "what": f"control = {name} over [{lo:.4g}, {hi:.4g}]", "spec": s})
    if ctrl and isinstance(ctrl.get("range"), list) and len(ctrl["range"]) == 2:
        lo, hi = float(ctrl["range"][0]), float(ctrl["range"][1])
        w = 10 * abs(hi - lo)
        for label, (a, b) in (("wider", (lo - w, hi + w)), ("above", (hi, hi + w)), ("below", (lo - w, lo))):
            s = dict(spec)
            s["control"] = {"name": ctrl["name"], "range": [a, b]}
            out.append({"move": "modify_protocol_order", "what": f"control {ctrl['name']} range {label}: [{a:.4g}, {b:.4g}]", "spec": s})
    for var, expr in drift.items():
        terms = split_terms(str(expr))
        if len(terms) < 2:
            continue
        for i, t in enumerate(terms):
            s = dict(spec)
            d = dict(drift)
            d[var] = join_terms(terms[:i] + terms[i + 1:])
            s["drift"] = d
            out.append({"move": "project_operator_apparatus", "what": f"drop '{t}' from d{var}/dt", "spec": s})
    if donor is not None:
        out.extend(_donor_terms(spec, donor))
    if not check_bounded:
        return out
    kept = []
    for c in out:
        ok, why = bounded(c["spec"])
        if ok:
            kept.append(c)
        elif rejected is not None:
            rejected.append({"move": c["move"], "what": c["what"], "why": why})
    return kept


def _donor_terms(spec: Dict, donor: Dict) -> List[Dict]:
    """Additive terms of the donor's equations, moved to the specification's variable in the same position, the
    donor's parameters set to one and a new parameter g_add as the coefficient."""
    drift = dict(spec.get("drift") or {})
    params = dict(spec.get("parameters") or {})
    ctrl = spec.get("control") if isinstance(spec.get("control"), dict) else None
    fvars = list((donor.get("drift") or {}).keys())
    cvars = list(drift.keys())
    out: List[Dict] = []
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
                continue  # a symbol that is neither a variable nor a known function remains
            var = mapping[fv]
            s = dict(spec)
            d = dict(drift)
            d[var] = f"{drift[var]} + g_add*({mapped})"
            s["drift"] = d
            pr = dict(params)
            pr["g_add"] = 1.0
            s["parameters"] = pr
            s["material"] = [p for p in pr if not (ctrl and ctrl.get("name") == p)]
            out.append({"move": "add_operator_apparatus",
                        "what": f"add g_add*({mapped}) to d{var}/dt (from {donor.get('name', 'the donor')})", "spec": s})
    return out


# ------------------------------------------------------------------------------------------------ judgement
def judge(spec: Dict, seed: int = 0, targets=TARGETS) -> Dict[str, Dict]:
    """The write targets derived on `spec` in this process, without law checks: status, obstruction and seconds."""
    from . import codiscovery as cd
    try:
        real = mspec.load(spec)
    except Exception as e:  # noqa: BLE001
        return {t: {"status": "not read", "obstruction": f"{type(e).__name__}: {str(e)[:100]}", "seconds": 0.0} for t in targets}
    out = {}
    for t in targets:
        t0 = time.time()
        try:
            row = cd.TARGETS[t](real, np.random.default_rng(seed), check_law=False)
            out[t] = {"status": str(row.get("status", "")), "obstruction": str(row.get("obstruction", ""))[:160],
                      "seconds": round(time.time() - t0, 1)}
        except Exception as e:  # noqa: BLE001
            out[t] = {"status": "not derived", "obstruction": f"{type(e).__name__}: {str(e)[:100]}",
                      "seconds": round(time.time() - t0, 1)}
    return out


def judge_isolated(spec: Dict, seed: int = 0, budget: int = 120, targets=TARGETS) -> Dict[str, Dict]:
    """`judge` with every target in a child process of its own under a wall-clock limit (an exception raised inside a
    solver's callback would end this interpreter; a child that runs over is stopped and reported)."""
    import tempfile
    from .cli import _child_step
    settings = {"seed": seed, "law": False, "trajectories": 0, "samples": 0, "trajectories_completion": 0}
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "candidate.json"
        path.write_text(json.dumps(spec), encoding="utf-8")
        for t in targets:
            child, why = _child_step(path, settings, t, budget)
            if child is not None and child.get("not_read"):
                out[t] = {"status": "not read", "obstruction": child["not_read"], "seconds": 0.0}
            elif child is not None and t in child.get("targets", {}):
                row = child["targets"][t]
                out[t] = {"status": str(row.get("status", "")), "obstruction": str(row.get("obstruction", row.get("note", "")))[:160],
                          "seconds": row.get("seconds", 0.0)}
            else:
                out[t] = {"status": "not derived", "obstruction": why, "seconds": float(budget)}
    return out


def summary_of(der: Dict) -> Tuple[set, Dict[str, str]]:
    reached = {t for t, r in der.items() if str(r.get("status", "")).startswith("reached")}
    classes = {t: obstruction_class(r.get("obstruction", ""))["kind"] for t, r in der.items()
               if not str(r.get("status", "")).startswith("reached")}
    return reached, classes


def outcome_key(der: Dict) -> Tuple[Tuple[str, ...], Tuple[Tuple[str, str], ...]]:
    reached, classes = summary_of(der)
    return tuple(sorted(reached)), tuple(sorted(classes.items()))


def closes(composed_key, known_keys: List) -> Tuple[bool, str]:
    """Whether a composition's outcome closes on the known outcomes (the specification's and the single steps')."""
    r_c, c_c = composed_key
    if composed_key in known_keys:
        return True, "equals a known outcome"
    known_reached = {t for r, _ in known_keys for t in r}
    new_targets = [t for t in r_c if t not in known_reached]
    if new_targets:
        return False, "reaches " + ", ".join(new_targets) + " that no single step reaches"
    known_classes = {(t, k) for _, cl in known_keys for t, k in cl}
    new_classes = [f"{t}: {k}" for t, k in c_c if (t, k) not in known_classes and k not in UNDECIDED]
    if new_classes:
        return False, "new obstruction class " + "; ".join(new_classes)
    if any(k in UNDECIDED for _, k in c_c):
        return True, "undecided (the derivation did not resolve a class)"
    return True, "a combination of known classes"


# ------------------------------------------------------------------------------------------------ the run
def run(spec: Dict, donor: Optional[Dict] = None, depth: int = 1, seed: int = 0, budget: int = 120,
        max_candidates: int = 60, singles: int = 8, seconds: float = 1500.0, isolate: bool = True,
        log=None) -> Dict[str, object]:
    """Single moves on `spec` (depth 1), compositions of two with their closure (depth 2) and a third move on the open
    ones (depth 3)."""
    def derive(s: Dict) -> Dict:
        return judge_isolated(s, seed, budget) if isolate else judge(s, seed)

    t0 = time.time()
    base = derive(spec)
    b_reached, b_classes = summary_of(base)
    base_key = outcome_key(base)
    rejected: List[Dict] = []
    cands = candidates(spec, donor, rejected)
    rep: Dict[str, object] = {"name": spec.get("name", ""), "field": spec.get("field", ""), "base": base,
                              "candidates": len(cands), "rejected": rejected, "steps": [], "judged": [],
                              "compositions": [], "depth": depth}
    single_keys: Dict[str, Tuple] = {}
    for c in cands[:max_candidates]:
        if time.time() - t0 > seconds:
            rep["stopped"] = f"the time limit of {seconds:.0f} s was reached after {len(rep['judged'])} candidates"
            break
        der = derive(c["spec"])
        single_keys[c["what"]] = outcome_key(der)
        reached, classes = summary_of(der)
        new = sorted(reached - b_reached)
        changed = {t: [b_classes.get(t), classes.get(t)] for t in classes if t in b_classes and classes[t] != b_classes[t]}
        entry = {"move": c["move"], "what": c["what"], "reaches": new, "obstruction_changed": changed, "derivation": der}
        rep["judged"].append({"move": c["move"], "what": c["what"], "reaches": new, "obstruction_changed": changed})
        if new or changed:
            rep["steps"].append(entry)
        if log:
            log(f"  {c['move']}: {c['what']} -> reaches {new or '-'}; changes {changed or '-'}")
    if depth < 2:
        rep["seconds"] = round(time.time() - t0, 1)
        return rep
    steps = {s["what"] for s in rep["steps"]}
    first = [c for c in cands if c["what"] in steps] + [c for c in cands if c["what"] not in steps]
    first = first[:max(singles, len(steps))]
    for c1 in first:
        if time.time() - t0 > seconds:
            rep["stopped"] = f"the time limit of {seconds:.0f} s was reached during the compositions"
            break
        key1 = single_keys.get(c1["what"]) or outcome_key(derive(c1["spec"]))
        single_keys[c1["what"]] = key1
        for c2 in candidates(c1["spec"], donor)[:12]:
            if time.time() - t0 > seconds:
                break
            if c2["what"] == c1["what"] or (c1["move"] == "add_protocol_order" and c2["move"] == "add_protocol_order"):
                continue
            if c1["move"] == "modify_protocol_order" and c2["move"] == "modify_protocol_order":
                continue  # a range widened twice runs away from the model's scales
            der2 = derive(c2["spec"])
            key2 = outcome_key(der2)
            known = [base_key, key1] + ([single_keys[c2["what"]]] if c2["what"] in single_keys else [])
            ok, why = closes(key2, known)
            entry = {"moves": [[c1["move"], c1["what"]], [c2["move"], c2["what"]]], "closes": ok, "why": why,
                     "outcome": {"reached": list(key2[0]), "classes": dict(key2[1])}, "derivation": der2}
            if log:
                log(f"  {c1['what']} THEN {c2['what']}: {'closes' if ok else 'OPEN'} ({why})")
            if not ok and depth >= 3:
                thirds = []
                for c3 in candidates(c2["spec"], donor)[:6]:
                    if time.time() - t0 > seconds:
                        break
                    if c3["what"] in (c1["what"], c2["what"]) or (
                            c3["move"] == "modify_protocol_order" and "modify_protocol_order" in (c1["move"], c2["move"])):
                        continue
                    key3 = outcome_key(derive(c3["spec"]))
                    ok3, why3 = closes(key3, known + [key2])
                    thirds.append({"move": [c3["move"], c3["what"]], "closes": ok3, "why": why3,
                                   "outcome": {"reached": list(key3[0]), "classes": dict(key3[1])}})
                entry["third_steps"] = thirds
            rep["compositions"].append(entry)
    rep["open"] = sum(1 for e in rep["compositions"] if not e["closes"])
    rep["seconds"] = round(time.time() - t0, 1)
    return rep


def markdown(rep: Dict[str, object]) -> str:
    base = rep["base"]
    lines = [f"# Typed moves on {rep['name']}", "",
             "Base: " + "; ".join(f"{t}: {r['status']}" + (f" ({r['obstruction'][:70]})" if r.get("obstruction") and not str(r["status"]).startswith("reached") else "")
                                for t, r in base.items()) + ".", "",
             f"{rep['candidates']} candidates kept, {len(rep['rejected'])} rejected as unbounded, {len(rep['judged'])} judged, "
             f"{len(rep['steps'])} steps (a move that reaches a target the base does not, or changes an obstruction class)."
             + (f" {rep['stopped']}." if rep.get("stopped") else ""), ""]
    if rep["steps"]:
        lines += ["| Move | Change | Reaches | Obstruction class changed |", "| --- | --- | --- | --- |"]
        for s in rep["steps"]:
            ch = "; ".join(f"{t}: {a} → {b}" for t, (a, b) in s["obstruction_changed"].items())
            lines.append(f"| {s['move']} | {s['what']} | {', '.join(s['reaches']) or ''} | {ch} |")
        lines.append("")
    if rep["rejected"]:
        lines += ["Rejected as unbounded:", ""] + [f"- {r['move']}: {r['what']} ({r['why']})" for r in rep["rejected"][:20]] + [""]
    if rep.get("compositions"):
        opens = [e for e in rep["compositions"] if not e["closes"]]
        lines += [f"## Compositions of two moves: {len(rep['compositions'])}, {len(opens)} open", ""]
        for e in opens[:40]:
            lines.append(f"- {e['moves'][0][1]} THEN {e['moves'][1][1]} → {e['why']}"
                         + (f"; third steps: " + "; ".join(f"{t['move'][1]} ({'closes' if t['closes'] else 'open'})" for t in e["third_steps"]) if e.get("third_steps") else ""))
        if not opens:
            lines.append("Every composition closes on the single steps.")
        lines.append("")
    return "\n".join(lines)
