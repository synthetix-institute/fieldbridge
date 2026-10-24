"""The record of a drift specification: every quantity the rules measure on it, in one schema.

A record is the unit the constructor judges and reports. It is assembled from six steps, each of which can run in a
process of its own under a time limit (the `memory record` command does that), so that a slow or failing step costs
its own entry and not the record:

  symmetric-write, threshold-write, phase-locking
      the three write targets as `memory codiscover` derives them: reached or obstructed, the obstruction, the
      derivation's word and class, the law constant with its error and the expected constant;
  write_points
      every point along the control where a stable state loses stability or appears, with the normal form there
      (kind, critical dimension, the coefficients of the reduced drift, the bias along the critical mode and the slope
      of the critical rate), the terms that obstruct a symmetric write and what a completion would keep; the order of
      the symmetry group of the drift and the number of stable states at the declared parameters;
  completion
      the completion rung at every write point (`memory complete`): the law of the completed model, its certificate,
      the response constant and offset in one critical direction, the write class after a hold at the end value of the
      control (linear-stage write, equilibrium write, nonlinear stage) and the fraction of states that change during
      that hold, which is the retention;
  attractor
      for a realization without a stable state at its declared parameters, the kind of motion (chaotic, cycle,
      quasi-periodic, drift, converges, unbounded) from the Lyapunov exponents and the recurrence of the motion
      (attractor.attractor), and the laws of chaotic loss and of blow-up where they apply (chaos_blowup).

The record carries the hash of the specification and of the implementation, so that two records are comparable only
when both agree. Nothing in a record establishes that the model describes a physical system or that a result is new.
"""
from __future__ import annotations

import hashlib
import json
import platform
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from . import spec as mspec

SCHEMA = "fieldbridge-record/1"
WRITE_TARGETS = ("symmetric-write", "threshold-write", "phase-locking")
STEPS = WRITE_TARGETS + ("write_points", "completion", "attractor")


def default_settings() -> Dict[str, object]:
    return {"seed": 0, "law": True, "trajectories": 400, "samples": 4000, "trajectories_completion": 2000}


# ------------------------------------------------------------------------------------------------ helpers
def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items() if k != "final"}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


def timed(fn, *args, **kwargs):
    """fn(*args, **kwargs) with the seconds it took and the exception, if any, as a note. No alarm interrupts it: an
    exception raised inside a solver's callback ends the interpreter, so time limits are wall-clock limits on the
    process that runs a step."""
    t0 = time.time()
    try:
        return fn(*args, **kwargs), time.time() - t0, ""
    except Exception as e:  # noqa: BLE001
        return None, time.time() - t0, f"{type(e).__name__}: {str(e)[:120]}"


def implementation_hash() -> str:
    here = Path(__file__).resolve().parent
    return hashlib.sha256(b"".join(p.read_bytes() for p in sorted(here.glob("*.py")))).hexdigest()


def versions() -> Dict[str, str]:
    out = {"python": platform.python_version(), "numpy": np.__version__}
    try:
        import scipy
        out["scipy"] = scipy.__version__
    except Exception:  # noqa: BLE001
        pass
    return out


def law_summary(law) -> Dict:
    if not isinstance(law, dict):
        return {"status": str(law)[:160]}
    out = {k: law[k] for k in ("status", "law_kind", "critical_dimension", "dynamical_critical_dimension", "second_rate",
                               "slow_modes", "certified", "sector_certified", "obstruction_in_nonlinear_stage", "write_class")
           if k in law}
    if isinstance(law.get("stable_states"), list):
        out["stable_states"] = len(law["stable_states"])
    rows = law.get("rows") or []
    if rows:   # the largest deviation of the two maps from the certificate, over the bias directions
        flow = [r.get("deviation_sigma") for r in rows if r.get("deviation_sigma") is not None]
        sect = [r.get("sector_deviation_sigma") for r in rows if r.get("sector_deviation_sigma") is not None]
        if flow:
            out["flow_map_deviation_sigma_max"] = float(max(flow))
        if sect:
            out["sector_deviation_sigma_max"] = float(max(sect))
    if isinstance(law.get("hold"), dict):
        out["hold"] = {k: law["hold"].get(k) for k in ("hold_time", "switching_fraction")}
    od = law.get("one_direction")
    if isinstance(od, dict):   # response constant and offset of a one-dimensional write, law and certificate
        out["one_direction"] = {src: [{k: r[k] for k in ("bias_floors", "response_constant", "offset") if k in r}
                                      for r in rows] for src, rows in od.items() if isinstance(rows, list)}
    return out


def one_direction_text(law: Dict) -> str:
    od = law.get("one_direction") or {}
    parts = []
    for src in ("law", "certificate"):
        rows = od.get(src) or []
        if rows:
            r = max(rows, key=lambda x: x.get("bias_floors", 0))
            parts.append(f"{src}: response {r['response_constant']:.2f}, offset {r['offset']:+.2f}")
    return "; ".join(parts)


def completion_summary(entry: Dict) -> Dict:
    out = {"completion": entry.get("completion"), "kind": entry.get("kind"), "write_point": entry.get("write_point")}
    if "law" in entry:
        out["law"] = law_summary(entry["law"])
    if "unfoldings" in entry:
        out["unfoldings"] = [{"second": u.get("second"), "success": u.get("success"),
                              **({"cusp": u["cusp"], "kind_along_sweep": u.get("kind_along_sweep")} if u.get("success") else {}),
                              **({"law": law_summary(u["law"])} if "law" in u else {})} for u in entry["unfoldings"]]
    return out


# ------------------------------------------------------------------------------------------------ steps
def _load(spec):
    return mspec.load(spec if isinstance(spec, dict) else json.loads(Path(spec).read_text(encoding="utf-8")))


def step(spec: Dict, name: str, settings: Optional[Dict] = None) -> Dict[str, object]:
    """One step of the record on a specification (a dict). The result carries `seconds` and, where the step could not
    be done, `note`."""
    settings = {**default_settings(), **(settings or {})}
    if name not in STEPS:
        raise ValueError(f"unknown step {name!r}; the steps are {STEPS}")
    rng = np.random.default_rng(settings["seed"])
    try:
        real = _load(spec)
    except Exception as e:  # noqa: BLE001
        return {"status": "not read", "note": f"{type(e).__name__}: {str(e)[:120]}", "seconds": 0.0}
    if name in WRITE_TARGETS:
        from . import codiscovery as cd
        res, secs, note = timed(cd.TARGETS[name], real, rng, check_law=settings["law"], n_traj=settings["trajectories"])
        if not isinstance(res, dict):
            return {"status": "not derived", "note": note or "no result", "seconds": round(secs, 1)}
        keep = {k: res[k] for k in ("status", "obstruction", "obstruction_short", "word", "class") if k in res}
        const = res.get("law_constant") or {}
        if const.get("constant") is not None:
            keep["law_constant"] = {k: const[k] for k in ("constant", "stderr", "expected", "offset", "offset_stderr", "offset_resolved") if k in const}
        keep["seconds"] = round(secs, 1)
        return jsonable(keep)
    if name == "write_points":
        res, secs, note = timed(_write_points, real, rng)
        if res is None:
            return {"status": "not derived", "note": note, "seconds": round(secs, 1)}
        res["seconds"] = round(secs, 1)
        return jsonable(res)
    if name == "completion":
        from . import completion as cm
        if not real.control:
            return {"status": "no control parameter", "write_points": [], "seconds": 0.0}
        res, secs, note = timed(cm.complete, real, rng, n_samples=settings["samples"],
                                n_traj=settings["trajectories_completion"], hold_always=True)
        if not isinstance(res, dict):
            return {"status": "not derived", "note": note or "no result", "write_points": [], "seconds": round(secs, 1)}
        return jsonable({"status": "derived", "write_points": [completion_summary(e) for e in res.get("write_points", [])],
                         "seconds": round(secs, 1)})
    if name == "attractor":
        res, secs, note = timed(_attractor, real, rng)
        if res is None:
            return {"status": "not derived", "note": note, "seconds": round(secs, 1)}
        res["seconds"] = round(secs, 1)
        return jsonable(res)
    raise AssertionError(name)


def _write_points(real, rng) -> Dict[str, object]:
    from . import analysis as an
    from . import completion as cm
    from . import construct, predict
    out: Dict[str, object] = {"status": "derived"}
    syms = predict.symmetries(real, np.random.default_rng(1))
    out["symmetry_order"] = 1 + sum(1 for g in syms if "permutation" in g)
    stable = cm.stable_states(real, np.random.default_rng(2), {})
    out["stable_states_at_declared_parameters"] = len(stable)
    if not real.control:
        out["control"] = None
        out["note"] = "no control parameter: write points cannot be searched"
        out["points"] = []
        return out
    out["control"] = {"name": real.control, "range": list(real.control_range)}
    events = an.locate_writes(real, rng, n_starts=24 if real.carrier.dim <= 4 else 12)
    points = []
    for ev in events:
        nf = an.normal_form(real, ev["q"], real.control, ev["v"])
        ob = construct.obstruction(nf)
        p = {"control_value": float(ev["v"]), "state": list(map(float, ev["q"])), "kind": nf.get("kind"),
             "critical_dimension": int(nf.get("critical_dimension", 1) or 1),
             "coefficients": {k: nf[k] for k in ("a0", "a1", "a2", "a3", "a4", "a5") if k in nf},
             "bias_along_mode": nf.get("bias_along_mode"), "kappa_slope": nf.get("kappa_slope"),
             "second_rate": nf.get("second_rate"), "spectrum": nf.get("spectrum"),
             "critical_block_rank": nf.get("critical_block_rank"),
             "dynamical_critical_dimension": an.dynamical_critical_dimension(nf, cm.RATE_DEFAULT, cm.SLOW_FACTOR),
             "obstruction": ob.get("terms", {}), "completion": cm.completion_kind(nf, ob)}
        if nf.get("leading_term"):
            p["leading_term"] = nf["leading_term"]
        points.append(p)
    out["points"] = points
    return out


def _attractor(real, rng) -> Dict[str, object]:
    from . import attractor as at
    from . import chaos_blowup as cb
    from . import completion as cm
    stable = cm.stable_states(real, np.random.default_rng(2), {})
    if stable:
        return {"status": "stable states", "stable_states": len(stable)}
    out: Dict[str, object] = {"status": "measured", "stable_states": 0}
    out["motion"] = at.attractor(real, rng)
    kind = out["motion"].get("kind")
    if kind == "converges":   # a stable state the random starts of the search did not reach
        out["status"] = "converges to a state outside the sampled box"
    if kind in ("chaotic", "unbounded", "drift"):   # a drift can be a slow escape: the blow-up law decides
        out["laws"] = cb.characterize(real, np.random.default_rng(3))
    return out


# ------------------------------------------------------------------------------------------------ assembly
def assemble(spec: Dict, parts: Dict[str, Dict], source: str = "") -> Dict[str, object]:
    """The record of `spec` from the results of its steps (missing steps are recorded as not run)."""
    carrier = spec.get("carrier") or {}
    control = spec.get("control") if isinstance(spec.get("control"), dict) else None
    rec: Dict[str, object] = {
        "schema": SCHEMA, "name": spec.get("name", ""), "field": spec.get("field", ""), "kind": spec.get("kind"),
        "source": source, "spec_sha256": mspec.input_hash(spec),
        "carrier": {"kind": carrier.get("kind"), "dimension": len(carrier.get("variables") or [])},
        "control": {"name": control["name"], "range": control.get("range")} if control else None,
        "material": list(spec.get("material") or []),
        "targets": {t: parts.get(t, {"status": "not run"}) for t in WRITE_TARGETS},
        "write_points": parts.get("write_points", {"status": "not run"}),
        "completion": parts.get("completion", {"status": "not run"}),
        "attractor": parts.get("attractor", {"status": "not run"}),
        "implementation_sha256": implementation_hash(), "versions": versions(),
    }
    if control is None:
        rec["note"] = "no control parameter declared: the write targets and the completion cannot be searched"
    wp = rec["write_points"]
    if isinstance(wp, dict):
        rec["symmetry_order"] = wp.get("symmetry_order")
        rec["stable_states_at_declared_parameters"] = wp.get("stable_states_at_declared_parameters")
    # retention: the fraction of states that change during a hold as long as the sweep, per completed write point
    ret = []
    for e in (rec["completion"].get("write_points") or []) if isinstance(rec["completion"], dict) else []:
        laws = ([e["law"]] if isinstance(e.get("law"), dict) else []) + [u["law"] for u in e.get("unfoldings", []) if isinstance(u.get("law"), dict)]
        for L in laws:
            if isinstance(L.get("hold"), dict) and L["hold"].get("switching_fraction") is not None:
                ret.append({"completion": e.get("completion"), "write_class": L.get("write_class"),
                            "switching_fraction": L["hold"]["switching_fraction"], "hold_time": L["hold"].get("hold_time")})
    rec["retention"] = ret
    rec["summary"] = summary_line(rec)
    return rec


def summary_line(rec: Dict) -> str:
    """One line a person can read: the writes reached, the completions and their classes, the motion."""
    t = rec["targets"]
    reached = [k for k in WRITE_TARGETS if str(t.get(k, {}).get("status", "")).startswith("reached")]
    parts = [("reaches " + ", ".join(reached)) if reached else "reaches no write target"]
    comp = rec.get("completion", {})
    classes = []
    for e in comp.get("write_points", []) if isinstance(comp, dict) else []:
        laws = ([e["law"]] if isinstance(e.get("law"), dict) else []) + [u["law"] for u in e.get("unfoldings", []) if isinstance(u.get("law"), dict)]
        for L in laws:
            if L.get("write_class"):
                classes.append(f"{e.get('completion')}: {L['write_class']}")
    if classes:
        parts.append("completions: " + "; ".join(classes))
    att = rec.get("attractor", {})
    if isinstance(att, dict) and att.get("status") == "measured":
        parts.append("motion without a stable state: " + str(att.get("motion", {}).get("kind")))
    return "; ".join(parts)


def record(spec, settings: Optional[Dict] = None, steps=STEPS, source: str = "") -> Dict[str, object]:
    """The whole record in this process (no time limits; `memory record` isolates the steps)."""
    spec_d = spec if isinstance(spec, dict) else json.loads(Path(spec).read_text(encoding="utf-8"))
    parts = {s: step(spec_d, s, settings) for s in steps}
    return assemble(spec_d, parts, source or (str(spec) if not isinstance(spec, dict) else ""))


# ------------------------------------------------------------------------------------------------ text
def _target_cell(t: Dict) -> str:
    st = str(t.get("status", ""))
    if st.startswith("reached"):
        const = t.get("law_constant") or {}
        if const.get("constant") is not None:
            return f"{st}; {const['constant']:.3f} ± {const.get('stderr', 0.0):.3f}"
        return st
    if st in ("not derived", "not run", "not read"):
        return f"{st} ({t.get('note', '')})" if t.get("note") else st
    return f"{st}: {t.get('obstruction_short') or str(t.get('obstruction', ''))[:70]}"


def markdown(records: List[Dict]) -> str:
    lines = ["# Records", "",
             "One row per specification: the write targets with their law constants, the write points along the control "
             "with their normal forms, the completions with their classes and retention, and the motion of a realization "
             "without a stable state. The schema is `fieldbridge-record/1`.", "",
             "| Specification | Symmetric write | Threshold write | Phase locking | Write points | Completions | Retention (switching) | Motion |",
             "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in records:
        t = r["targets"]
        wp = r.get("write_points", {})
        pts = wp.get("points") if isinstance(wp, dict) else None
        wtext = (", ".join(f"{p['kind'].split(':')[0]} at {p['control_value']:.3g}" for p in pts) if pts
                 else (wp.get("note") or wp.get("status", "") if isinstance(wp, dict) else ""))
        comp = r.get("completion", {})
        ctext = ""
        if isinstance(comp, dict):
            ctext = "; ".join(f"{e.get('completion')}: {(e.get('law') or {}).get('write_class') or (e.get('law') or {}).get('status', '')}"
                              for e in (comp.get("write_points") or []))
            if not ctext:   # no write point, or a step that was stopped or could not run
                ctext = comp.get("status", "") + (f" ({comp['note']})" if comp.get("note") else "")
        ret = r.get("retention") or []
        rtext = ", ".join(f"{x['switching_fraction']:.2f}" for x in ret)
        att = r.get("attractor", {})
        mtext = (att.get("motion", {}).get("kind", "") if att.get("status") == "measured" else
                 f"{att.get('stable_states', '')} stable states" if att.get("status") == "stable states" else att.get("status", ""))
        lines.append(f"| {r['name'][:40]} | {_target_cell(t['symmetric-write'])} | {_target_cell(t['threshold-write'])} | "
                     f"{_target_cell(t['phase-locking'])} | {wtext[:80]} | {ctext[:90]} | {rtext} | {mtext} |")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------ the predictions table
def load_records(path) -> List[Dict]:
    """Records from a `records.jsonl` or `records.json` file."""
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    if p.suffix == ".jsonl":
        return [json.loads(line) for line in text.splitlines() if line.strip()]
    data = json.loads(text)
    return data["records"] if isinstance(data, dict) else data


def target_cell(t: Dict) -> str:
    """The cell of a write target in the table of predictions: its status with the law constant, or what obstructs it."""
    st = str(t.get("status", ""))
    if st.startswith("reached"):
        const = t.get("law_constant") or {}
        if const.get("constant") is not None:
            return (f"{st}; law constant {const['constant']:.3f} ± {const.get('stderr', 0.0):.3f}"
                    + (f" (expected {const['expected']:.3f})" if const.get("expected") is not None else ""))
        return st
    if st in ("not derived", "not run", "not read"):
        return f"{st} ({t.get('note', '')})" if t.get("note") else st
    return f"{st}: {t.get('obstruction_short') or str(t.get('obstruction', ''))[:90]}"


def class_text(L: Dict) -> str:
    """The class a failed certificate leaves: the equilibrium write (the choice exchanged with the noise during a hold)
    or the nonlinear stage."""
    wc = L.get("write_class")
    if wc == "equilibrium write":
        h = L.get("hold") or {}
        return (f" (equilibrium write: {100 * h.get('switching_fraction', 0):.0f} % change state during a hold of "
                f"{h.get('hold_time', 0):.3g})")
    if wc == "nonlinear stage" or (L.get("obstruction_in_nonlinear_stage") and wc is None):
        return " (nonlinear stage)"
    return ""


def completion_cell(c: Dict) -> str:
    """The completions of a record or of a predictions row, with their laws, certificates, classes and one-direction
    constants."""
    if not c:
        return ""
    if c.get("note") and not c.get("write_points"):
        return f"not run ({c['note']})"
    parts = []
    for e in c.get("write_points", []):
        kind = e.get("completion", "")
        if "law" in e:
            L = e["law"]
            cert = "" if "certified" not in L else (", certified" if L["certified"] else ", not certified")
            od = one_direction_text(L)
            parts.append(f"{kind}: {L.get('status', '')}{cert}{class_text(L)}" + (f"; {od}" if od else ""))
        elif "unfoldings" in e:
            done = [u for u in e["unfoldings"] if u.get("success")]
            if done:
                u = done[0]
                L = u.get("law", {})
                cert = "" if "certified" not in L else (", certified" if L["certified"] else ", not certified")
                od = one_direction_text(L)
                parts.append(f"{kind}: cusp unfolded along {u['second']} ({u.get('kind_along_sweep', '')} along the "
                             f"sweep); law {L.get('status', '')}{cert}{class_text(L)}" + (f"; {od}" if od else ""))
            else:
                parts.append(f"{kind}: no material parameter cancels the asymmetry")
        else:
            parts.append(kind)
    if parts:
        return "; ".join(parts)
    return c.get("status", "") if c.get("status") not in (None, "derived") else "no write point along the control"


def predictions_rows(records: List[Dict]) -> List[Dict]:
    """The rows of the table of the write targets (`memory predictions` on drift specifications) read off the records:
    the same keys, so that the table is a view of the records and not a second derivation."""
    rows = []
    for r in records:
        src = str(r.get("source") or "")
        rows.append({"spec": Path(src).name if src else r.get("name", ""), "name": r.get("name", ""),
                     "field": r.get("field", ""), "kind": r.get("kind"), "rule": "write point", "source": src,
                     "targets": r.get("targets", {}), "completion": r.get("completion", {}),
                     "spec_sha256": r.get("spec_sha256"), "summary": r.get("summary", "")})
    return rows


def _laws_of(e: Dict) -> List[Dict]:
    return ([e["law"]] if isinstance(e.get("law"), dict) else []) + [u["law"] for u in e.get("unfoldings", []) if isinstance(u.get("law"), dict)]


def predictions_counts(records: List[Dict]) -> Dict[str, object]:
    """What the table reports in its first paragraph: targets reached, obstructed symmetric writes, certified
    completions, write classes by specification, steps stopped at a time limit."""
    reached = {t: 0 for t in WRITE_TARGETS}
    obstructed_symmetric, certified, stopped = 0, 0, 0
    classes: Dict[str, List[str]] = {}
    for r in records:
        t = r.get("targets", {})
        for k in WRITE_TARGETS:
            st = str(t.get(k, {}).get("status", ""))
            reached[k] += st.startswith("reached")
            stopped += "time limit" in str(t.get(k, {}).get("note", ""))
        obstructed_symmetric += str(t.get("symmetric-write", {}).get("status", "")).startswith("obstructed")
        for step in ("write_points", "completion", "attractor"):
            stopped += "time limit" in str((r.get(step) or {}).get("note", ""))
        comp = r.get("completion", {})
        for e in (comp.get("write_points") or []) if isinstance(comp, dict) else []:
            laws = _laws_of(e)
            certified += any(L.get("certified") for L in laws)
            for L in laws:
                if L.get("write_class"):
                    classes.setdefault(L["write_class"], []).append(Path(str(r.get("source") or r.get("name", ""))).stem)
    return {"specifications": len(records), "reached": reached, "obstructed_symmetric": obstructed_symmetric,
            "completions_certified": certified, "stopped": stopped,
            "write_classes": {k: sorted(set(v)) for k, v in sorted(classes.items())}}


def predictions_markdown(records: List[Dict], title: str = "The write rules on the recorded specifications",
                         intro: str = "") -> str:
    """The table of the write targets and completions, derived from records."""
    c = predictions_counts(records)
    n = c["specifications"]
    classes = "; ".join(f"{k}: {', '.join(v)}" for k, v in c["write_classes"].items())
    head = (f"Derived from the records (schema `{SCHEMA}`) of {n} specification(s). {c['reached']['symmetric-write']} reach "
            f"the symmetric write along their declared control, {c['reached']['threshold-write']} the threshold write and "
            f"{c['reached']['phase-locking']} phase locking; {c['obstructed_symmetric']} have an obstructed symmetric write; "
            f"the completion rung certified a law at {c['completions_certified']} write point(s); {c['stopped']} step(s) "
            f"stopped at a time limit." + (f" Write classes of the completed models ({classes})." if classes else ""))
    lines = [f"# {title}", "", (intro + " " if intro else "") + head, "",
             "## Drift equations and networks: the write targets", "",
             "Each target is read from the record of the specification, where `memory record` derived it as `memory codiscover` "
             "does: reached (the model's write point carries the mechanism, with the law constant where the law was checked), or "
             "obstructed, with what obstructs it. The completion column lists the completion rung at every write point of the "
             "record (`memory complete`, with the hold that fixes the write class): the law of the completed model, its "
             "certificate from the full simulation, and the response constant and offset in one critical direction.", "",
             "| Specification | Field | Symmetric write | Threshold write | Phase locking | Completion at the write points |",
             "| --- | --- | --- | --- | --- | --- |"]
    for r in predictions_rows(records):
        t = r["targets"]
        lines.append(f"| {r['spec']} | {r['field']} | {target_cell(t.get('symmetric-write', {}))} | "
                     f"{target_cell(t.get('threshold-write', {}))} | {target_cell(t.get('phase-locking', {}))} | "
                     f"{completion_cell(r.get('completion', {}))} |")
    return "\n".join(lines) + "\n"
