"""``fieldbridge memory ...``: build, predict and transfer memory from a specification file.

  predict  SPEC                structure only: what the realization can store and how (no simulation)
  card     SPEC                the calculated memory card, the structural predictions and their comparison
  attach   --from SPEC --to SPEC   which properties of a memory are preserved, changed or obstructed on another carrier
  design   SPEC --free PARAM   the setting of PARAM that restores a symmetric (pitchfork) write
  phase    SPEC                memory in the phase of an oscillating realization (phase response, retention, writing)
  regimes                      information about a write for kappa > 0, = 0 and < 0
  field    SPEC                memory in a field: conservation and dimension set the law of loss
  dynamics SPEC                chaotic loss of the initial state and blow-up in finite time, with their laws
  hysterons SPEC               return to a turning point of a slow drive: signs of couplings and drive against subloops
  counterpart SPEC             the balanced counterpart of a hysteron network and the states measured against it
  gallery                      every example: cards, figures and an HTML gallery (index.html)
  loops                        loops of genes, spins and rotors: holonomy against simulation
  networks                     random networks of genes, spins and rotors: structure against simulation
  changes  SPEC                changes of a model and pairs of changes, judged by the symmetric and threshold writes
  codiscover [SPEC ...]        one mechanism (--target symmetric-write, threshold-write or phase-locking) derived in
                               realizations from different fields: derivations, obstructions, invariants (default:
                               every example)

Contributing a material (CONTRIBUTING.md):
  new NAME                     a template specification that loads and runs, with placeholders to replace
  check SPEC [SPEC ...]        required and recommended checks, the structural predictions and a quick card
  record PATH                  the record of a specification or a directory: targets, write points, completions, motion
  changes SPEC                 changes of the specification (a term removed or added from a donor, another control, a wider range), judged by the write targets; pairs and their closure
  catalog                      the table of materials in docs/materials.md, calculated from the specifications

Every command writes <name>.json and <name>.md to --out-dir with the question, the input hash, the hashes of the
implementation, library versions, and the evidence boundary (novelty is never established by the program).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
from typing import Dict

BOUNDARY = ("Calculated and structural statements about a supplied model. They do not establish that the model "
            "describes a physical system, and they do not establish novelty; both need comparison with experiment "
            "and with the literature.")


def _provenance(spec: Dict = None) -> Dict[str, object]:
    here = Path(__file__).resolve().parent
    impl = hashlib.sha256(b"".join(p.read_bytes() for p in sorted(here.glob("*.py")))).hexdigest()
    versions = {"python": platform.python_version()}
    for mod in ("numpy", "scipy", "sympy", "matplotlib"):
        try:
            versions[mod] = __import__(mod).__version__
        except ImportError:
            versions[mod] = None
    out = {"implementation_sha256": impl, "versions": versions, "novelty_established": False,
           "source_alignment_verified": False, "evidence_boundary": BOUNDARY}
    if spec is not None:
        from .spec import input_hash
        out.update(question=spec.get("question"), input_sha256=input_hash(spec), assumptions=spec.get("assumptions"),
                   provenance=spec.get("provenance"))
    return out


def _write(out_dir: Path, name: str, report: Dict, markdown: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{name}.json").write_text(json.dumps(report, indent=2, default=float) + "\n", encoding="utf-8")
    (out_dir / f"{name}.md").write_text(markdown, encoding="utf-8")


# Display names for the keys of the reports (the JSON keys are kept for compatibility).
VERDICT_NAMES = {"stores": "Stable states", "writes": "Write point", "constructs": "Writing mechanism",
                 "holds": "Retention law", "lock": "Writing protocols"}
PREDICTION_NAMES = {"write": "Writing", "holding": "Retention", "lock": "Retention and rewriting",
                    "oscillation": "Oscillation", "multistability": "Multistability", "network": "Loops (holonomy)"}


SIGNATURE_NAMES = {"reciprocal": "Reciprocal couplings", "write": "Writing", "holding": "Retention",
                   "oscillation": "Oscillation", "multistability": "Multistability",
                   "control_role": "Role of the control", "continuous_symmetry": "Continuous symmetry (zero mode)",
                   "satisfiable": "All loops satisfiable", "reflect_to_align": "Ratio of reflection to alignment coupling",
                   "carrier": "Carrier"}


def _value(v) -> str:
    """A signature entry in words: booleans as yes/no, an undecided entry as such."""
    if v is None:
        return "not decided by the structure"
    if isinstance(v, bool):
        return "yes" if v else "no"
    return str(v)


def _md_predictions(pred: Dict) -> str:
    lines = []
    for key, p in pred.get("predictions", {}).items():
        body = p["prediction"] if isinstance(p["prediction"], str) else json.dumps(p["prediction"], default=float)
        lines.append(f"- **{PREDICTION_NAMES.get(key, key)}**: {body}  \n  rule: {p['rule']}; uses: {p['uses']}; "
                     f"falsified by: {p['falsifier']}")
    return "\n".join(lines)


def cmd_predict(args) -> int:
    import numpy as np
    from . import predict, spec
    real = spec.load(args.spec)
    pred = predict.predict(real, np.random.default_rng(args.seed))
    report = {"command": "predict", **_provenance(real.spec), "structure": pred}
    md = (f"# Predictions from structure\n\n{real.spec['question']}\n\n" + _md_predictions(pred)
          + f"\n\n{BOUNDARY}\n")
    _write(Path(args.out_dir), "predict", report, md)
    print(json.dumps({"out_dir": args.out_dir, "predictions": sorted(pred.get("predictions", {}))}, indent=2))
    return 0


def cmd_card(args) -> int:
    import numpy as np
    from . import discovery, predict, spec
    real = spec.load(args.spec)
    pred = predict.predict(real, np.random.default_rng(args.seed))
    card = discovery.evaluate(real, np.random.default_rng(args.seed + 1), quick=args.quick)
    card["id"] = Path(args.spec).stem
    comp = predict.compare(pred, card)
    report = {"command": "card", **_provenance(real.spec), "structure": pred, "card": card, "comparison": comp}
    v = card["verdict"]
    md = (f"# Memory card: {card['name']}\n\n{real.spec['question']}\n\n## Calculated\n\n"
          + "\n".join(f"- **{VERDICT_NAMES.get(k, k)}**: {val}" for k, val in v.items())
          + "\n\n## Predicted from structure\n\n" + _md_predictions(pred)
          + "\n\n## Agreement\n\n" + "\n".join(f"- {k}: {c}" for k, c in comp.items())
          + f"\n\n{BOUNDARY}\n")
    out = Path(args.out_dir)
    _write(out, "card", report, md)
    if not args.no_figure:
        from .visual import card_figure
        card_figure(card, out / "card.png")
    print(json.dumps({"out_dir": args.out_dir, "verdict": v, "all_consistent": comp.get("all_consistent")}, indent=2))
    return 0


def _signature(pred: Dict) -> Dict[str, object]:
    P = pred.get("predictions", {})
    net = P.get("network", {}).get("prediction", {}) if "network" in P else {}
    return {"reciprocal": pred.get("reciprocal"),
            "write": P.get("write", {}).get("prediction"),
            "holding": P.get("holding", {}).get("prediction"),
            "oscillation": P.get("oscillation", {}).get("prediction"),
            "multistability": P.get("multistability", {}).get("prediction"),
            "control_role": pred.get("control_role"),
            "continuous_symmetry": bool(pred.get("continuous_symmetries")) or bool(net.get("continuous_symmetry")),
            "satisfiable": net.get("satisfiable")}


def _slots(real) -> Dict[str, object]:
    """What the realization is made of, as opposed to what its memory does: carrier, transports, coupling ratio."""
    import numpy as np
    out = {"carrier": real.carrier.describe()}
    net = getattr(real, "network", None)
    if net is not None and net["type"] == "rotor":
        kinds = set()
        ratios = []
        for e in net["edges"]:
            terms = {e[2]: 1.0} if isinstance(e[2], str) else e[2]
            kinds |= set(terms)
            if terms.get("align") and terms.get("reflect"):
                ratios.append(terms["reflect"] / terms["align"])
        out["transports"] = sorted(kinds)
        if ratios:
            out["reflect_to_align"] = f"median {float(np.median(ratios)):.3g} (range {min(ratios):.3g}-{max(ratios):.3g})"
    return out


def cmd_attach(args) -> int:
    import numpy as np
    from . import predict, spec
    src, tgt = spec.load(args.source), spec.load(args.target)
    ps = predict.predict(src, np.random.default_rng(args.seed))
    pt = predict.predict(tgt, np.random.default_rng(args.seed))
    ss, st = _signature(ps), _signature(pt)
    kept = {k: ss[k] for k in ss if ss[k] == st[k]}
    changed = {k: {"source": ss[k], "target": st[k]} for k in ss if ss[k] != st[k]}
    a, b = _slots(src), _slots(tgt)
    materials = {k: {"source": a.get(k), "target": b.get(k)} for k in set(a) | set(b) if a.get(k) != b.get(k)}
    moves = [
        {"move": "detach", "acts_on": "Omega", "result": "the memory signature of the source, independent of its carrier"},
        {"move": "attach", "acts_on": "Xi", "result": f"the target carrier: {tgt.carrier.describe()}"},
        {"move": "close", "acts_on": "C", "result": tgt.closure or "target closure not stated"},
        {"move": "observe", "acts_on": "R", "result": tgt.observable or "target observable not stated"},
        {"move": "execute", "acts_on": "P", "result": "write protocol of the target (field, pulse or sweep of its control)"},
        {"move": "falsify", "acts_on": "I_op", "result": "each changed entry names a prediction that a simulation "
                                                         "or an experiment on the target must confirm"},
    ]
    report = {"command": "attach", **_provenance(tgt.spec), "source": src.spec.get("question"),
              "source_signature": ss, "target_signature": st, "kept": kept, "changed": changed,
              "carrier_and_material": materials, "moves": moves}
    md = (f"# Attaching a memory\n\nFrom: {src.spec['question']}\n\nTo: {tgt.spec['question']}\n\n"
          "## Memory signature kept\n\n"
          + "\n".join(f"- {SIGNATURE_NAMES.get(k, k)}: {_value(v)}" for k, v in kept.items())
          + "\n\n## Memory signature changed\n\n"
          + ("\n".join(f"- {SIGNATURE_NAMES.get(k, k)}: {_value(v['source'])} -> {_value(v['target'])}"
                       for k, v in changed.items()) or "- nothing")
          + "\n\n## Carrier and material changed\n\n"
          + ("\n".join(f"- {SIGNATURE_NAMES.get(k, k)}: {_value(v['source'])} -> {_value(v['target'])}"
                       for k, v in materials.items()) or "- nothing")
          + "\n\n## Moves\n\n" + "\n".join(f"- {m['move']} ({m['acts_on']}): {m['result']}" for m in moves)
          + f"\n\n{BOUNDARY}\n")
    _write(Path(args.out_dir), "attach", report, md)
    print(json.dumps({"out_dir": args.out_dir, "kept": sorted(kept), "changed": sorted(changed),
                      "carrier_and_material": sorted(materials)}, indent=2))
    return 0


def _md_design(i: int, e: Dict, control: str, free: str) -> str:
    wp = e["write_point"]
    state = ", ".join(f"{x:.4g}" for x in wp["state"])
    ob = "; ".join(f"{k} = {v:.3g}" for k, v in e["obstruction"].items()) or "none"
    lines = [f"## Write point {i}: {control} = {wp[control]:.4g}, state ({state})", "",
             f"- kind: {e['kind']}", f"- obstruction to the supercritical pitchfork: {ob}"]
    c = e.get("cancellation")
    if c is None:
        lines.append("- no cancellation needed")
    elif not c.get("success"):
        lines.append(f"- no setting of {free} removes the obstruction within the search range"
                     + (f": {c['bound_note']}" if c.get("bound_note") else ""))
    else:
        d = c["direction"]
        cstate = ", ".join(f"{x:.4g}" for x in c["state"])
        lines += [f"- the obstruction vanishes at {free} = {c['value']:.4g}, {control} = {c['control']:.4g}, "
                  f"state ({cstate})",
                  f"- sweep through this point along ({control}, {free}) proportional to "
                  f"({d[control]:.3g}, {d[free]:.3g})",
                  f"- kind along the sweep: {c['kind_along_sweep']}"]
    for s in (c or {}).get("bounds_reached") or ():
        lines.append(f"- the fold followed toward the bound {s['text']} stops there, at {free} = {s[free]:.4g}, "
                     f"{control} = {s[control]:.4g}")
    return "\n".join(lines)


def cmd_design(args) -> int:
    import numpy as np
    from . import analysis as an, construct, spec
    real = spec.load(args.spec)
    if args.free not in real.params:
        raise SystemExit(f"--free {args.free!r} is not a parameter of the specification")
    if not real.control:
        raise SystemExit("the specification declares no control parameter to cross a write point")
    rng = np.random.default_rng(args.seed)
    events = an.locate_writes(real, rng, n_starts=16)
    results = []
    for ev in events:
        nf = an.normal_form(real, ev["q"], real.control, ev["v"])
        ob = construct.obstruction(nf)
        entry = {"write_point": {real.control: ev["v"], "state": ev["q"]}, "kind": nf["kind"], "obstruction": ob["terms"]}
        if ob["asymmetric"] or ob["biased"]:
            c = construct.cancel_asymmetry(real, nf, args.free)
            entry["cancellation"] = {k: c.get(k) for k in ("success", "value", "control", "state", "direction")}
            entry["cancellation"].update({k: c[k] for k in ("bounds_reached", "cusp_beyond_bound", "bound_note")
                                          if k in c})
            if c.get("success"):
                entry["cancellation"]["kind_along_sweep"] = c["normal_form_along_sweep"]["kind"]
        results.append(entry)
    report = {"command": "design", **_provenance(real.spec), "free_parameter": args.free, "write_points": results}
    md = (f"# Parameter setting that restores a symmetric write\n\n{real.spec['question']}\n\n"
          f"Free parameter: {args.free}. Control parameter: {real.control}.\n\n"
          + "\n\n".join(_md_design(i + 1, e, real.control, args.free) for i, e in enumerate(results))
          + f"\n\n{BOUNDARY}\n")
    _write(Path(args.out_dir), "design", report, md)
    print(json.dumps({"out_dir": args.out_dir, "write_points": len(results)}, indent=2))
    return 0


COMPLETION_NAMES = {"none": "no obstruction (the law is derived as a control)",
                    "critical space": "keep every critical direction",
                    "leading order": "keep the leading term of higher order",
                    "second parameter": "keep a second parameter and sweep through the cusp",
                    "distant branch": "keep the distant branch to which the state jumps",
                    "phase": "keep the phase (phase-locking target)"}


def _md_complete(spec_question: str, rep: Dict) -> str:
    lines = ["# The law at an obstruction", "", spec_question, ""]
    for i, e in enumerate(rep["write_points"], 1):
        wp = e["write_point"]
        lines += [f"## Write point {i}", "", f"- kind: {e['kind']}",
                  f"- obstruction: " + ("; ".join(f"{k} = {v:.3g}" for k, v in e["obstruction"].items()) or "none"),
                  f"- completion: {COMPLETION_NAMES.get(e['completion'], e['completion'])}"]
        laws = [("", e["law"])] if "law" in e else [(f" (second parameter {u['second']})", u.get("law"))
                                                      for u in e.get("unfoldings", []) if u.get("law")]
        for tag, L in laws:
            lines.append(f"- law{tag}: {L.get('status')}")
            if L.get("status") != "derived":
                continue
            lines.append(f"  - critical directions {L['critical_dimension']}, noise {L['noise']:.3g} ({L['noise_rule']}), "
                         f"rate {L['rate']}")
            for src, rows in (L.get("one_direction") or {}).items():
                for x in rows:
                    lines.append(f"  - {src}: bias {x['bias_floors']} floors, response constant "
                                 f"{x['response_constant']:.3f} (pi^1/4 = 1.331 for the pitchfork), offset {x['offset']:.3f}")
            if L["critical_dimension"] > 1:
                lines.append(f"  - stable states at the end of the sweep: {len(L['stable_states'])}; probabilities by bias "
                             "direction in the critical space (law / certificate):")
                for r in L["rows"]:
                    lines.append(f"    - {[round(d, 2) for d in r['direction']]} ({r['bias_floors']} floors): "
                                 f"{[round(p, 3) for p in r['law']]} / {[round(p, 3) for p in r.get('certificate', [])]}")
            if "certified" in L:
                lines.append("  - certificate: " + ("the law agrees with the full stochastic dynamics" if L["certified"]
                                                     else L.get("obstruction_in_nonlinear_stage", "not certified")))
        lines.append("")
    return "\n".join(lines) + f"\n{BOUNDARY}\n"


def cmd_audit(args) -> int:
    """The specific ways the symmetric-write law of a specification could be wrong, each settled by its own test."""
    import numpy as np
    from . import audit, spec
    real = spec.load(args.spec)
    res = audit.audit_symmetric_write(real, np.random.default_rng(args.seed), n_traj=args.trajectories,
                                      replicates=args.replicates, checks=args.checks)
    report = {"command": "audit", **_provenance(real.spec), "question": real.spec.get("question"),
              "settings": {"trajectories": args.trajectories, "replicates": args.replicates, "checks": args.checks},
              **json.loads(json.dumps(res, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))}
    _write(Path(args.out_dir), "audit", report, audit.markdown(real.spec.get("name", ""), res) + f"\n{BOUNDARY}\n")
    print(json.dumps({"out_dir": args.out_dir, "verdict": res["verdict"],
                      "checks": {x["check"]: x["status"] for x in res.get("ledger", [])}}, indent=2))
    return 0


def cmd_complete(args) -> int:
    import numpy as np
    from . import completion, spec
    real = spec.load(args.spec)
    rep = completion.complete(real, np.random.default_rng(args.seed), rate=args.rate, noise=args.noise,
                              floors=tuple(args.floors), n_samples=args.samples, n_traj=args.trajectories,
                              certify=not args.no_certificate, gamma=args.gamma, curvature=args.curvature)

    def clean(o):
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items() if k != "final"}
        if isinstance(o, list):
            return [clean(v) for v in o]
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        return o
    rep = clean(rep)
    report = {"command": "complete", **_provenance(real.spec), "question": real.spec.get("question"),
              "settings": {"rate": args.rate, "gamma": args.gamma, "noise": args.noise, "floors": args.floors,
                           "samples": args.samples, "trajectories": args.trajectories,
                           "curvature": args.curvature}, **rep}
    _write(Path(args.out_dir), "complete", report, _md_complete(real.spec.get("question", ""), rep))
    print(json.dumps({"out_dir": args.out_dir, "write_points": [
        {"completion": e["completion"],
         "certified": [L.get("certified") for L in ([e["law"]] if "law" in e else
                                                     [u.get("law", {}) for u in e.get("unfoldings", [])])]}
        for e in rep["write_points"]]}, indent=2))
    return 0


def cmd_changes(args) -> int:
    import sys
    from . import changes, spec
    real = spec.load(args.spec)
    donor = json.loads(Path(args.donor).read_text(encoding="utf-8")) if args.donor else None
    res = changes.changes(real.spec, budget=args.budget, seed=args.seed, pairs=args.pairs, first=args.first,
                          kinds=tuple(args.kinds), only=args.only, donor=donor, isolate=not args.no_isolate,
                          progress=(lambda m: print(m, file=sys.stderr, flush=True)) if args.verbose else None)
    out = Path(args.out_dir)
    if args.save:
        d = out / "changed"
        d.mkdir(parents=True, exist_ok=True)
        rows = [r for r in res["singles"] if r["step"]] + [p for p in res.get("pairs", []) if not p["closes"]]
        for i, r in enumerate(rows):
            (d / f"change_{i:03d}.json").write_text(json.dumps(r["spec"], indent=2) + "\n", encoding="utf-8")
    for r in res["singles"]:
        r.pop("spec", None)
    for p in res.get("pairs", []):
        p.pop("spec", None)
    report = {"command": "changes", **_provenance(real.spec), "question": real.spec.get("question"),
              "settings": {"budget": args.budget, "pairs": args.pairs, "first": args.first, "kinds": args.kinds,
                           "only": args.only, "donor": args.donor, "isolated": not args.no_isolate}, **res}
    _write(out, "changes", report, changes.markdown(real.spec.get("name", ""), res) + f"\n{BOUNDARY}\n")
    print(json.dumps({"out_dir": args.out_dir, "steps": res["steps"],
                      "reached": sorted({t for r in res["singles"] for t in r["reaches"]}),
                      "pairs": len(res.get("pairs", [])), "open": res.get("open")}, indent=2))
    return 0


def cmd_loops(args) -> int:
    import numpy as np
    from . import compose
    rows = compose.run_loops(np.random.default_rng(args.seed))
    ok = sum(r["consistent"] for r in rows)
    report = {"command": "loops", **_provenance(), "loops": rows, "consistent": ok, "total": len(rows)}
    _write(Path(args.out_dir), "loops", report, f"# Loops across fields\n\n{ok} of {len(rows)} as predicted.\n\n{BOUNDARY}\n")
    if not args.no_figure:
        from .visual import loops_figure
        loops_figure(rows, Path(args.out_dir) / "loops.png")
    print(json.dumps({"out_dir": args.out_dir, "consistent": ok, "total": len(rows)}, indent=2))
    return 0


def cmd_networks(args) -> int:
    import numpy as np
    from . import networks
    rows = networks.benchmark(np.random.default_rng(args.seed), {"gene": args.count, "spin": args.count,
                                                                 "rotor": args.count})
    by = {}
    for r in rows:
        by.setdefault(r["kind"], [0, 0])
        by[r["kind"]][0] += r["consistent"]
        by[r["kind"]][1] += 1
    report = {"command": "networks", **_provenance(), "networks": rows,
              "summary": {k: {"consistent": a, "total": b} for k, (a, b) in by.items()}}
    md = "# Random networks: structure against simulation\n\n" + "\n".join(
        f"- {k}: {a} of {b} as predicted" for k, (a, b) in by.items()) + f"\n\n{BOUNDARY}\n"
    _write(Path(args.out_dir), "networks", report, md)
    print(json.dumps(report["summary"], indent=2))
    return 0


def cmd_phase(args) -> int:
    import numpy as np
    from . import phase, spec
    real = spec.load(args.spec)
    rng = np.random.default_rng(args.seed)
    try:
        cyc = phase.limit_cycle(real, rng)
    except ValueError as error:
        report = {"command": "phase", **_provenance(real.spec), "refused": str(error)}
        _write(Path(args.out_dir), "phase", report, f"# Memory in the phase\n\nRefused: {error}\n\n{BOUNDARY}\n")
        print(json.dumps({"out_dir": args.out_dir, "refused": str(error)}, indent=2))
        return 1
    kick = [args.kick] + [0.0] * (real.carrier.dim - 1)
    prc = phase.response_curve(real, cyc, rng, kick=kick, n=24)
    hold = phase.hold(real, cyc, rng, shift=args.shift, D=args.noise, n_traj=200, n_periods=40)
    lock = phase.lock(real, cyc, rng, strengths=[0.01, 0.02, 0.05, 0.1, 0.2], target=args.shift, D=args.noise)
    report = {"command": "phase", **_provenance(real.spec), "cycle": cyc, "response_curve": prc, "hold": hold,
              "lock": lock}
    md = (f"# Memory in the phase\n\n{real.spec['question']}\n\n- period {cyc['period']:.4g}\n"
          f"- a kick of {args.kick:g} on the first coordinate shifts the phase by {min(prc['shift']):.3g} to "
          f"{max(prc['shift']):.3g} rad, depending on when it arrives\n"
          f"- a written shift of {args.shift:.3g} rad stays at {hold['separation'][-1]:.3g} rad after "
          f"{hold['time'][-1]:.4g} time units; the phase variance grows as 2 D_phi t with D_phi = "
          f"{hold['phase_diffusion']:.3g} (Law 2)\n"
          f"- ratio of retention time to writing time against the write strength: log-log slope "
          f"{lock['log_log_slope']:.3f} (linear, as expected along a flat direction)\n\n"
          f"{BOUNDARY}\n")
    _write(Path(args.out_dir), "phase", report, md)
    print(json.dumps({"out_dir": args.out_dir, "period": cyc["period"], "phase_diffusion": hold["phase_diffusion"],
                      "lock_slope": lock["log_log_slope"]}, indent=2))
    return 0


def cmd_regimes(args) -> int:
    import numpy as np
    from . import regimes
    r = regimes.three_regimes(np.random.default_rng(args.seed), kappa=args.kappa, s0=args.s0, D=args.noise, n=args.n)
    report = {"command": "regimes", **_provenance(),
              "input": {"kappa": args.kappa, "s0": args.s0, "noise": args.noise, "n": args.n, "seed": args.seed}, **r}
    rows = r["regimes"]
    t = np.array(r["times"])
    pick = [0, len(t) // 2, len(t) - 1]
    lines = [f"- {k} (kappa = {v['kappa']:g}): information " + ", ".join(
        f"{v['gaussian_information'][i]:.3g} at t = {t[i]:.3g}" for i in pick) for k, v in rows.items()]
    md = ("# Information about a write in the three regimes of kappa\n\n"
          "I(t) = (1/2) ln[1 + (kappa s0^2/D) / (exp(2 kappa t) - 1)]\n\n" + "\n".join(lines)
          + f"\n\nAt kappa < 0 the information stops falling at (1/2) ln(1 + W) = {r['plateau_gaussian']:.3g}, "
            f"W = |kappa| s0^2 / D = {r['W']:.3g}: the expansion outruns the noise. The sign of the observable is "
            f"kept with probability Phi(W^1/2) = {r['binary_plateau']['accuracy']:.3f} (simulated "
            f"{rows['expanding']['simulation']['accuracy'][-1]:.3f} at t = {t[-1]:.3g}).\n\n{BOUNDARY}\n")
    _write(Path(args.out_dir), "regimes", report, md)
    print(json.dumps({"out_dir": args.out_dir, "W": r["W"], "plateau": r["plateau_gaussian"]}, indent=2))
    return 0


def cmd_noise_ratio(args) -> int:
    import numpy as np
    from . import noise_ratio as nr, spec
    if args.law or not args.spec:
        md = nr.law_markdown() + f"\n{BOUNDARY}\n"
        report = {"command": "noise-ratio", **_provenance(), "law": {
            "gammas": list(nr.GAMMAS), "continuous": list(nr.CONTINUOUS), "first_order": list(nr.FIRST_ORDER),
            "early_gammas": list(nr.EARLY_GAMMAS), "early": list(nr.EARLY)}}
        _write(Path(args.out_dir), "noise_ratio", report, md)
        print(md)
        return 0
    real = spec.load(args.spec)
    if not real.control:
        raise SystemExit("the specification declares no control parameter to sweep")
    rng = np.random.default_rng(args.seed)
    nf = nr.symmetric_point(real, rng)
    if nf is None:
        report = {"command": "noise-ratio", **_provenance(real.spec),
                  "refused": "no symmetric critical point (continuous or first-order pitchfork) along the control; "
                             "a fold may be unfolded to a cusp first (memory design)"}
        _write(Path(args.out_dir), "noise_ratio", report, f"# Write law across the noise ratio\n\n{report['refused']}\n")
        print(json.dumps({"out_dir": args.out_dir, "refused": report["refused"]}, indent=2))
        return 1
    gammas = [float(g) for g in args.gammas.split(",")]
    rep = nr.sweeps(real, nf, rng, gammas=gammas, q=args.q, n_traj=args.trajectories)
    report = {"command": "noise-ratio", **_provenance(real.spec), **rep}
    _write(Path(args.out_dir), "noise_ratio", report, nr.markdown(rep) + f"\n{BOUNDARY}\n")
    print(json.dumps({"out_dir": args.out_dir, "first_order": rep["first_order"],
                      "constants": {f"{r['gamma']:g}": round(r["constant"], 4) for r in rep["rows"]}}, indent=2))
    return 0


def cmd_dynamics(args) -> int:
    """memory dynamics: chaotic loss of the initial state and blow-up in finite time, each with its law
    (fieldbridge.memory.chaos_blowup)."""
    import numpy as np
    from . import chaos_blowup as cb
    from .spec import load
    spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
    real = load(spec)
    rng = np.random.default_rng(args.seed)
    ch = cb.chaos(real, rng, starts=args.starts)
    # the blow-up test integrates each start until it escapes; it runs when a start escaped in the transient
    bu = (cb.blowup(real, rng, starts=args.starts) if ch.get("escaped")
          else {"verdict": "no start escapes during the transient of the chaos test"})
    report = {"command": "dynamics", **_provenance(spec), "chaos": ch, "blowup": bu}
    lines = [f"# Chaos and blow-up of {spec.get('name', 'the model')}", "", str(spec.get("question", "")), "",
             f"- chaos: {ch.get('verdict')}"]
    pr = ch.get("predictability") or {}
    if pr.get("law") is not None:
        lines.append(f"  - predictability: separation times {['%.3g' % t if t is not None else '-' for t in pr['times']]} "
                     f"at delta = {pr['deltas']}; slope x lambda_1 = {pr['law']:.3f} (the law gives 1)")
    lines.append(f"- blow-up: {bu.get('verdict')}")
    if bu.get("balances"):
        lines.append("  - dominant balances (growth exponents p): "
                     + ", ".join(f"{b['growth']:.3g}" for b in bu["balances"]))
    sc = bu.get("escape_time_scaling")
    if sc:
        lines.append(f"  - escape time of {'weighted' if sc['weighted'] else 'uniformly'} scaled starts: slope {sc['slope']:.3f} (expected "
                     f"{sc['expected']:.3f})")
    _write(Path(args.out_dir), "dynamics", report, "\n".join(lines) + f"\n\n{BOUNDARY}\n")
    print(json.dumps({"out_dir": args.out_dir, "chaos": ch.get("verdict"), "blowup": bu.get("verdict")}, indent=2))
    return 0


def cmd_field(args) -> int:
    import numpy as np
    from . import fields
    spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
    model = fields.from_spec(spec)
    r = fields.run(model, np.random.default_rng(args.seed), pairs=args.pairs)
    report = {"command": "field", **_provenance(spec), **r}
    pred = r["prediction"]
    if pred["law"] == "power":
        ren = r["renormalized"]
        body = (f"- predicted from conservation, dimension and the shape of the write: SNR ~ t^{pred['exponent']:g}"
                f"\n- exact linear field, infinite-lattice limit ({r['exponent_check']['L']} sites per side): exponent "
                f"{r['exponent_check']['exponent']:.3f}"
                f"\n- nonlinear field (g = {model.g:g}) against the linear curve with noise and clock renormalized by the "
                f"site variance (r = T/<phi^2> = {ren['r']:.3f}, no fitted parameter): rms log residual "
                f"{ren['rms_log_residual']:.3f} for t >= {ren['t_min']:g}"
                f"\n- {pred['global']}; information kept by the total of {model.L ** model.d} sites: "
                f"{r['total_charge_information']:.3g} nats")
    else:
        h = r["hartree"]
        fmt = lambda v: "n/a" if v is None else f"{v:.3g}"  # noqa: E731
        body = (f"- predicted: exponential loss at rate {pred['rate']:.3g} ({pred['rate_is']}; no conservation)"
                f"\n- rates over t = {fmt(r['window'][0])} to {fmt(r['window'][1])}: linear field "
                f"{fmt(r['rate_linear'])}; simulated field "
                f"{fmt(r['rate_simulated'])}; self-consistent mass m^2 = {h['mass2']:.3g} gives {fmt(h['rate'])}")
    md = f"# Memory in a field\n\n{spec['question']}\n\n{body}\n\n{BOUNDARY}\n"
    _write(Path(args.out_dir), "field", report, md)
    summary = {"out_dir": args.out_dir, "prediction": pred}
    if pred["law"] == "power":
        summary.update(exponent_exact=r["exponent_check"]["exponent"],
                       renormalized={k: r["renormalized"][k] for k in ("r", "rms_log_residual", "t_min")})
    else:
        summary.update(window=r["window"], rate_linear=r["rate_linear"], rate_simulated=r["rate_simulated"],
                       hartree=r["hartree"])
    print(json.dumps(summary, indent=2))
    return 0


def cmd_hysterons(args) -> int:
    import numpy as np
    from . import hysterons as hy
    spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
    if spec.get("kind") == "equations":
        return _hysterons_from_equations(args, spec)
    model = hy.from_spec(spec)
    c = hy.card(model, np.random.default_rng(args.seed), count=args.subloops)
    r = c["check"]
    pred = r["prediction"]
    if not args.no_figure:
        from .visual import hysteron_figure
        Path(args.out_dir).mkdir(parents=True, exist_ok=True)
        hysteron_figure(c, Path(args.out_dir) / "hysterons.png")
    report = {"command": "hysterons", **_provenance(spec), "elements": model.n, "class": c["class"], **r}
    names = {"least": "least stable first", "random": "random order", "least_nested": "least stable first, nested",
             "random_nested": "random order, nested"}
    rows = "\n".join(f"| {names[k]} | {v['subloops']} | {v['failed']} | {v['did_not_end']} | {v['max_differ']} |"
                     for k, v in r["subloops"].items())
    plaq = pred["frustrated_plaquettes"]
    md = (f"# Return to a turning point\n\n{spec.get('question', '')}\n\n"
          f"- structure: {pred['statement']}\n"
          f"- loops of couplings alone: {'none frustrated' if pred['couplings_balanced'] else 'some frustrated'}"
          + (f" (frustrated plaquettes: {plaq:.3f})" if plaq is not None else "")
          + f"; couplings whose loop through the drive is frustrated: {pred['bonds_frustrated_with_drive']} of "
          f"{pred['bonds']}\n- measured on {model.n} elements: {r['verdict']}\n\n"
          "| Relaxation of avalanches | Subloops | Did not return | Avalanche did not end | Most elements unlike at the "
          "return |\n| --- | --- | --- | --- | --- |\n" + rows + f"\n\n{BOUNDARY}\n")
    _write(Path(args.out_dir), "hysterons", report, md)
    print(json.dumps({"out_dir": args.out_dir, "guaranteed": pred["guaranteed"], "failed": r["failed"],
                      "total": r["total"], "verdict": r["verdict"]}, indent=2))
    return 0


def _hysterons_from_equations(args, spec: Dict) -> int:
    """memory hysterons on a material given by equations: its transition graph under the slow drive of the block
    slow_drive, the exact return test on it, and the reduction to hysterons (fieldbridge.memory.slow_drive)."""
    from . import slow_drive as sd
    from .spec import load
    real = load(spec)
    r = sd.check(spec, real, seed=args.seed)
    report = {"command": "hysterons", "from": "equations", **_provenance(spec), **r}
    _write(Path(args.out_dir), "hysterons", report, sd.markdown(spec, r) + f"\n{BOUNDARY}\n")
    print(json.dumps({"out_dir": args.out_dir, "states": r["count"], "return": r["return"],
                      "guaranteed": r["prediction"]["guaranteed"], "verdict": r["verdict"]}, indent=2))
    return 0


def cmd_counterpart(args) -> int:
    import numpy as np
    from . import counterpart as cp
    from . import hysterons as hy
    spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
    model = hy.from_spec(spec)
    r = cp.construct(model, np.random.default_rng(args.seed), max_states=args.max_states,
                     exact_up_to=args.exact_up_to, return_test=not args.no_return_test)
    report = {"command": "counterpart", **_provenance(spec), "elements": model.n, **r}

    def row(name, m):
        if m.get("censored"):
            return f"| {name} | more than {m['max_states']} | | | | |"
        if m.get("cycle"):
            return f"| {name} | an avalanche did not end | | | | |"
        ret = m.get("return")
        test = f"{ret['failed']} of {ret['tested']}" if ret else "not tested"
        return f"| {name} | {m['count']} | {m['C']:.2f} | {m['deep']:.2f} | {m['max_reversals']} | {test} |"

    fr = r["frustration"]
    md = (f"# Balanced counterpart\n\n{spec.get('question', '')}\n\n"
          f"- rule: {r['rule']}\n"
          f"- frustration index of the network with the drive: {fr['F']} "
          f"({'exact' if fr['exact'] else 'upper bound from local search'}); edges: {fr['edges_couplings']} "
          f"couplings, {fr['edges_drive']} actions of the drive\n"
          f"- reversed in the counterpart: couplings {fr['unsatisfied_couplings']}, drive on elements "
          f"{fr['unsatisfied_drive']}\n\n{r['statement']}\n\n"
          "| Material | States | Bits | Deep memory (bits) | Most reversals needed | Return test: failed excursions |\n"
          "| --- | --- | --- | --- | --- | --- |\n"
          + "\n".join([row("network", r["network"]), row("ferromagnetic reference (unique)", r["ferromagnetic_reference"]),
                       row("balanced counterpart (one minimizer)", r["balanced"]), row("uncoupled elements", r["uncoupled"])])
          + f"\n\n{BOUNDARY}\n")
    _write(Path(args.out_dir), "counterpart", report, md)
    summary = {"out_dir": args.out_dir, "F": fr["F"], "F_exact": fr["exact"]}
    summary.update({k: r[k] for k in ("gain", "loss", "typical_in_ensemble") if k in r})
    print(json.dumps(summary, indent=2))
    return 0


def cmd_order(args) -> int:
    """The order of writes measured against the uncoupled reference: the holonomy rule as a measurement."""
    from . import order as od
    spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
    net = od.from_spec(spec)
    r = od.construct(net, exact_up_to=args.exact_up_to)
    report = {"command": "order", **_provenance(spec), "elements": net.n, "writes": [w.name for w in net.writes], **r}
    net_m, ref_m = r["network"], r["uncoupled_reference"]
    rows = "\n".join(f"| {o} | {s or 'not certified'} | {ref_m['states'].get(o) or 'not certified'} |" for o, s in net_m["states"].items())
    md = (f"# Order of writes\n\n{spec.get('question', '')}\n\n- rule: {r['rule']}\n"
          f"- commutation index H = {r['index']['H']} ({'exact' if r['index']['exact'] else 'bound'}); "
          f"assignment that makes the writes commute: {r['index']['assignment']}\n"
          f"- density class {r['density_class']}, couplings {r['couplings']}; typical record in the census: "
          f"{r['typical_record_in_census']}\n\n{r['statement']}\n\n"
          "| Order | Final state of the network | Final state of the uncoupled reference |\n| --- | --- | --- |\n"
          + rows + f"\n\n{BOUNDARY}\n")
    _write(Path(args.out_dir), "order", report, md)
    print(json.dumps({"out_dir": args.out_dir, "H": r["index"]["H"], "distinct": net_m["distinct"], "split": net_m["split"],
                      "reference_distinct": ref_m["distinct"], "gain": r["gain"]}, indent=2))
    return 0


from .record import STEPS as STEPS_FOR_PARSER, WRITE_TARGETS, completion_summary as _completion_summary, jsonable as _jsonable  # noqa: E402
from .record import law_summary as _law_summary, one_direction_text as _one_direction_text, timed as _timed  # noqa: E402
from .record import target_cell as _drift_cell, completion_cell as _completion_cell  # noqa: E402


def _predict_drift(path: str, settings: Dict, progress: "Path | None" = None, only: "str | None" = None) -> Dict:
    """The write rules on one specification of kind equations or network: the three targets, and, where the symmetric
    write is obstructed, the completion rung (the law of the model that keeps what the obstruction names). With `only`,
    one write target or "completion" alone (the completion then runs whatever the symmetric write's status: the parent
    process that asks for it has decided). The row so far is written to `progress` after every step."""
    import numpy as np
    from . import codiscovery as cd
    from . import spec as mspec
    p = Path(path)
    spec = json.loads(p.read_text(encoding="utf-8"))
    row = {"spec": p.name, "name": spec.get("name", p.stem), "field": spec.get("field", ""), "kind": spec.get("kind"),
           "rule": "write point", "source": spec.get("provenance", {}).get("source", ""), "targets": {}}

    def save():
        if progress is not None:
            progress.write_text(json.dumps(_jsonable(row)), encoding="utf-8")

    try:
        real = mspec.load(spec)
    except Exception as e:  # noqa: BLE001
        row["not_read"] = f"{type(e).__name__}: {str(e)[:120]}"
        save()
        return row
    for target in (WRITE_TARGETS if only is None else [t for t in WRITE_TARGETS if t == only]):
        res, secs, note = _timed(cd.TARGETS[target], real, np.random.default_rng(settings["seed"]),
                                 check_law=settings["law"], n_traj=settings["trajectories"])
        if not isinstance(res, dict):
            row["targets"][target] = {"status": "not derived", "note": note or "no result", "seconds": round(secs, 1)}
        else:
            keep = {k: res[k] for k in ("status", "obstruction", "obstruction_short", "word", "class") if k in res}
            const = res.get("law_constant") or {}
            if const.get("constant") is not None:
                keep["law_constant"] = {k: const[k] for k in ("constant", "stderr", "expected", "offset", "offset_stderr", "offset_resolved") if k in const}
            keep["seconds"] = round(secs, 1)
            row["targets"][target] = keep
        save()
    sw = row["targets"].get("symmetric-write", {})
    if only == "completion" or (only is None and settings["complete"]
                                and str(sw.get("status", "")).startswith("obstructed")):
        from . import completion as cm
        row["completion"] = {"status": "running"}
        save()
        res, secs, note = _timed(cm.complete, real, np.random.default_rng(settings["seed"]),
                                 n_samples=settings["samples"], n_traj=settings["trajectories_completion"])
        row["completion"] = {"seconds": round(secs, 1), **({"note": note} if note else {}),
                             "write_points": [_completion_summary(e) for e in (res or {}).get("write_points", [])]}
        save()
    return _jsonable(row)


def _child_step(path: Path, settings: Dict, step: str, limit: int) -> "tuple[Dict | None, str]":
    """One step (a write target or the completion) of one specification in a child process with a wall-clock limit.
    Returns (the child's row or None, why it has none)."""
    import shutil
    import subprocess
    import sys
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp) / "spec"
        d.mkdir()
        shutil.copy(path, d / path.name)
        out = Path(tmp) / "out"
        cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "predictions", str(d), "--out-dir", str(out),
               "--no-isolate", "--only-step", step, "--seed", str(settings["seed"]),
               "--trajectories", str(settings["trajectories"]), "--samples", str(settings["samples"]),
               "--trajectories-completion", str(settings["trajectories_completion"])]
        if not settings["law"]:
            cmd.append("--no-law")
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=limit if limit > 0 else None)
        except subprocess.TimeoutExpired:
            return None, f"stopped at the time limit of {limit} s"
        except OSError as e:  # the interpreter could not be started
            return None, f"the derivation process could not start: {e}"
        report = out / "predictions.json"
        if res.returncode == 0 and report.exists():
            rows = json.loads(report.read_text(encoding="utf-8")).get("rows", [])
            if rows:
                return rows[0], ""
        tail = [ln for ln in (res.stderr or "").strip().splitlines() if ln.strip()]
        return None, (f"the derivation process ended with code {res.returncode}" + (f": {tail[-1][:160]}" if tail else ""))


def _predict_drift_isolated(path: str, settings: Dict) -> Dict:
    """_predict_drift with every step in a child process of its own under a wall-clock limit: `budget` seconds for
    each write target, `completion_budget` for the completion rung. A step that runs over is stopped and reported as
    such; a crash inside a numerical routine costs that step and not the run."""
    p = Path(path)
    spec = json.loads(p.read_text(encoding="utf-8"))
    row = {"spec": p.name, "name": spec.get("name", p.stem), "field": spec.get("field", ""), "kind": spec.get("kind"),
           "rule": "write point", "source": spec.get("provenance", {}).get("source", ""), "targets": {}}
    for target in WRITE_TARGETS:
        child, why = _child_step(p, settings, target, settings["budget"])
        if child is not None and child.get("not_read"):
            row["not_read"] = child["not_read"]
            return row
        if child is not None and target in child.get("targets", {}):
            row["targets"][target] = child["targets"][target]
        else:
            row["targets"][target] = {"status": "not derived", "note": why or "no result"}
    sw = row["targets"].get("symmetric-write", {})
    if settings["complete"] and str(sw.get("status", "")).startswith("obstructed"):
        child, why = _child_step(p, settings, "completion", settings["completion_budget"])
        if child is not None and "completion" in child:
            row["completion"] = child["completion"]
        else:
            row["completion"] = {"note": why or "no result", "write_points": []}
    if not any("seconds" in t for t in row["targets"].values()):
        row["failed"] = next((t.get("note", "") for t in row["targets"].values() if t.get("note")), "no step was derived")
    return row


def _record_child(path: Path, settings: Dict, step: str, limit: int) -> "tuple[Dict | None, str]":
    """One step of the record on one specification in a child process with a wall-clock limit (an exception raised
    inside a solver's callback would end this interpreter; a child that runs over is stopped and reported)."""
    import shutil
    import subprocess
    import sys
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp) / "spec"
        d.mkdir()
        shutil.copy(path, d / path.name)
        out = Path(tmp) / "out"
        cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "record", str(d), "--out-dir", str(out),
               "--no-isolate", "--only-step", step, "--seed", str(settings["seed"]),
               "--trajectories", str(settings["trajectories"]), "--samples", str(settings["samples"]),
               "--trajectories-completion", str(settings["trajectories_completion"])]
        if not settings["law"]:
            cmd.append("--no-law")
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=limit if limit > 0 else None)
        except subprocess.TimeoutExpired:
            return None, f"stopped at the time limit of {limit} s"
        except OSError as e:
            return None, f"the process could not start: {e}"
        report = out / "records.json"
        if res.returncode == 0 and report.exists():
            rows = json.loads(report.read_text(encoding="utf-8")).get("records", [])
            if rows and step in rows[0].get("steps", {}):
                return rows[0]["steps"][step], ""
        tail = [ln for ln in (res.stderr or "").strip().splitlines() if ln.strip()]
        return None, (f"the process ended with code {res.returncode}" + (f": {tail[-1][:160]}" if tail else ""))


def _record_one(path: str, settings: Dict, limits: Dict[str, int], isolate: bool, only_step: "str | None",
                steps=None) -> Dict:
    from . import record as rc
    p = Path(path)
    spec = json.loads(p.read_text(encoding="utf-8"))
    steps = [only_step] if only_step else list(steps or rc.STEPS)
    parts: Dict[str, Dict] = {}
    for st in steps:
        if isolate:
            res, why = _record_child(p, settings, st, limits.get(st, 0))
            parts[st] = res if res is not None else {"status": "not derived", "note": why}
        else:
            parts[st] = rc.step(spec, st, settings)
    if only_step:   # a child: the step alone, for its parent to assemble
        return {"spec": p.name, "steps": parts}
    return rc.assemble(spec, parts, source=str(p))


def cmd_record(args) -> int:
    """The record of every drift specification of a file or a directory (schema fieldbridge-record/1): the write
    targets, the write points with their normal forms, the completions with their classes and retention, and the
    motion of a realization without a stable state; every step in a child process under a time limit."""
    from . import record as rc
    out = Path(args.out_dir)
    if args.from_records:   # no derivation: the tables of an existing records file
        records = rc.load_records(args.from_records)
        stem = Path(args.from_records).stem
        out.mkdir(parents=True, exist_ok=True)
        table = rc.predictions_markdown(records, intro=args.intro or "") + f"\n{BOUNDARY}\n"
        (out / f"{stem}.md").write_text(table, encoding="utf-8")
        report = {"command": "record", "derived_from": str(args.from_records), **_provenance(),
                  "counts": rc.predictions_counts(records), "rows": rc.predictions_rows(records)}
        (out / f"{stem}.json").write_text(json.dumps(report, indent=1, default=float) + "\n", encoding="utf-8")
        print(json.dumps({"out_dir": args.out_dir, "derived_from": str(args.from_records), **rc.predictions_counts(records)}, indent=2))
        return 0
    if args.path is None:
        raise SystemExit("a specification file or a directory is needed (or --from-records FILE)")
    root = Path(args.path)
    if root.is_dir():
        paths = sorted(root.rglob("*.json") if args.recursive else root.glob("*.json"))
    else:
        paths = [root]
    paths = [p for p in paths if isinstance(json.loads(p.read_text(encoding="utf-8")), dict)
             and "schema" in json.loads(p.read_text(encoding="utf-8"))
             and json.loads(p.read_text(encoding="utf-8")).get("kind") in ("equations", "network")]
    settings = {"seed": args.seed, "law": not args.no_law, "trajectories": args.trajectories, "samples": args.samples,
                "trajectories_completion": args.trajectories_completion}
    limits = {t: args.budget for t in rc.WRITE_TARGETS}
    limits.update({"write_points": args.budget, "completion": args.completion_budget, "attractor": args.attractor_budget})
    isolate = not args.no_isolate
    steps = tuple(args.steps.split(",")) if args.steps else rc.STEPS
    unknown = [s for s in steps if s not in rc.STEPS]
    if unknown:
        raise SystemExit(f"unknown step(s) {unknown}; the steps are {list(rc.STEPS)}")
    if isolate and args.workers > 1:
        import concurrent.futures as cf
        with cf.ThreadPoolExecutor(max_workers=args.workers) as pool:
            records = list(pool.map(lambda p: _record_one(str(p), settings, limits, True, None, steps), paths))
    else:
        records = [_record_one(str(p), settings, limits, isolate, args.only_step, steps) for p in paths]
    out.mkdir(parents=True, exist_ok=True)
    report = {"command": "record", **_provenance(), "settings": {**settings, **limits, "isolated": isolate}, "records": records}
    (out / "records.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    with (out / "records.jsonl").open("w", encoding="utf-8") as fh:
        for r in records:
            fh.write(json.dumps(r) + "\n")
    if not args.only_step:
        (out / "records.md").write_text(rc.markdown(records) + f"\n{BOUNDARY}\n", encoding="utf-8")
        # the table of the write targets and completions, a view of the records
        table = rc.predictions_markdown(records) + f"\n{BOUNDARY}\n"
        (out / "predictions.md").write_text(table, encoding="utf-8")
        report = {"command": "record", "derived_from": "records.jsonl", **_provenance(),
                  "counts": rc.predictions_counts(records), "rows": rc.predictions_rows(records)}
        (out / "predictions.json").write_text(json.dumps(report, indent=1, default=float) + "\n", encoding="utf-8")
        print(json.dumps({"out_dir": args.out_dir, "records": len(records),
                          "summaries": [r.get("summary", "") for r in records][:12]}, indent=2))
    else:
        print(json.dumps({"out_dir": args.out_dir, "step": args.only_step, "records": len(records)}))
    return 0


def cmd_predictions(args) -> int:
    """What the rules predict for every specification in a directory. Hysteron networks: rule R8 (the frustration
    index with the drive, balance, the gain over the ferromagnetic reference and the return test where the
    enumeration is feasible, the typical gain of the sequence-memory ensemble for comparison). Bistable networks of
    kind order: the holonomy index of the writes. Quantum realizations: the closure rule. Drift equations and
    networks: the three write targets, each in a child process under a time limit, and the completion rung where
    the symmetric write is obstructed."""
    import numpy as np
    from . import counterpart as cp
    from . import hysterons as hy
    rows, lines, drift_paths, skipped = [], [], [], 0
    root = Path(args.directory)
    for path in sorted(root.rglob("*.json") if args.recursive else root.glob("*.json")):
        try:
            spec = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(spec, dict) or "schema" not in spec:
                continue   # not a specification (a survey or an index kept beside the models)
            if spec.get("kind") == "order":
                from . import order as od
                net = od.from_spec(spec)
                idx = od.commutation_index(net, exact_up_to=args.exact_up_to)
                typ = od.CENSUS.get("attractive" if all(s > 0 for s in net.signs) else "antiferromagnetic", {}).get(od.density_class(net))
                row = {"spec": path.name, "name": net.name, "elements": net.n, "rule": "holonomy", "H": idx["H"], "H_exact": idx["exact"],
                       "commuting": idx["H"] == 0, "prediction": "the writes commute: no order is recorded" if idx["H"] == 0 else
                       "the writes fail to commute: a record of the order is possible (necessary, not sufficient)",
                       "typical_record_in_census": typ, "source": spec.get("provenance", {}).get("source", "")}
                measured = ""
                if net.n <= args.enumerate_up_to:
                    r = od.construct(net, exact_up_to=args.exact_up_to)
                    row.update(distinct=r["network"]["distinct"], split=r["network"]["split"],
                               reference_distinct=r["uncoupled_reference"]["distinct"], gain=r["gain"])
                    measured = (f"distinct {r['network']['distinct']}, split {r['network']['split']}, reference "
                                f"{r['uncoupled_reference']['distinct']}, gain {r['gain']:.2f} bits")
                rows.append(row)
                lines.append(f"| {path.name} | {net.n} | order: H = {idx['H']} | {'yes' if idx['H'] == 0 else 'no'} (commute) | "
                             f"{'' if typ is None else f'{100 * typ:.0f} % record'} | {measured} |")
                continue
            if spec.get("schema") == "fieldbridge-quantum/1":
                from ..quantum import closure as qc
                from ..quantum import language as ql
                r = qc.construct(ql.load(spec))
                e = r.get("excess") or {}
                row = {"spec": path.name, "name": r["name"], "elements": r["sector_dimension"], "rule": "closure",
                       "closure_dimension": r["closure"]["dimension"], "closure_frequencies": len(r["closure"]["frequencies"]),
                       "index": r["index"], "rotating": r["index"] == 0, "rule_holds": r["rule_holds"],
                       "prediction": "the observable rotates with one frequency" if r["index"] == 0 else
                                     "the signal has several frequencies, all in the closure spectrum; removing "
                                     f"{r['index']} term(s) leaves one rotation" if r["index"] else r.get("obstruction", ""),
                       "relative_spread": e.get("relative_spread"), "mean_shift": e.get("mean_shift"),
                       "source": spec.get("provenance", {}).get("reference", "")}
                rows.append(row)
                measured = (f"{e['lines']} lines, relative spread {e['relative_spread']:.3f}, mean shift {e['mean_shift']:+.3g}"
                            if e else r.get("obstruction", ""))
                lines.append(f"| {path.name} | {r['sector_dimension']} | closure: {r['closure']['dimension']} operators, "
                             f"{len(r['closure']['frequencies'])} frequencies | index {'' if r['index'] is None else r['index']} | "
                             f"{'yes' if r['index'] == 0 else 'no'} (rotates) | {'holds' if r['rule_holds'] else 'violated'} | {measured} |")
            if spec.get("kind") in ("equations", "network") and spec.get("schema") == "fieldbridge-memory/1":
                drift_paths.append(path)
                continue
            if spec.get("kind") != "hysterons":
                skipped += 1
                continue
            model = hy.from_spec(spec)
        except Exception as e:  # noqa: BLE001
            lines.append(f"| {path.name} | | | | | | not read: {type(e).__name__} |")
            continue
        fr = cp.frustration(model, np.random.default_rng(args.seed), exact_up_to=args.exact_up_to)
        degree = 2.0 * len(model.undirected()) / max(model.n, 1)
        typ = cp.expected_gain(fr["F"], model.n, degree)
        row = {"spec": path.name, "name": model.name, "elements": model.n, "mean_degree": degree, "F": fr["F"],
               "F_exact": fr["exact"], "balanced": fr["F"] == 0,
               "prediction": "return to every turning point is exact" if fr["F"] == 0 else
                             "states beyond return-point memory are possible; the rule is sufficient, not necessary",
               "typical_gain_in_ensemble": typ["gain"], "source": spec.get("provenance", {}).get("source", "")}
        measured = ""
        if model.n <= args.enumerate_up_to:
            r = cp.construct(model, np.random.default_rng(args.seed), max_states=args.max_states, exact_up_to=args.exact_up_to)
            if "gain" in r:
                row.update(gain=r["gain"], loss=r["loss"], deep_gain=r["deep_gain"],
                           return_failed=r["network"]["return"]["failed"], return_tested=r["network"]["return"]["tested"],
                           reference_return_failed=r["ferromagnetic_reference"]["return"]["failed"])
                measured = (f"gain {r['gain']:.2f} bits, deep {r['deep_gain']:.2f}, return failures "
                            f"{r['network']['return']['failed']}/{r['network']['return']['tested']} "
                            f"(reference {r['ferromagnetic_reference']['return']['failed']})")
            else:
                measured = "enumeration stopped"
        rows.append(row)
        lines.append(f"| {path.name} | {model.n} | {degree:.1f} | {fr['F']}{'' if fr['exact'] else ' (bound)'} | "
                     f"{'yes' if fr['F'] == 0 else 'no'} | {typ['gain']:.2f} | {measured} |")
    # drift equations and networks: the write targets and the completion rung, one process per worker
    settings = {"seed": args.seed, "budget": args.budget, "law": not args.no_law, "trajectories": args.trajectories,
                "complete": not args.no_complete, "completion_budget": args.completion_budget, "samples": args.samples,
                "trajectories_completion": args.trajectories_completion}
    drift_rows = []
    if drift_paths:
        if args.no_isolate:
            progress = Path(args.progress_file) if args.progress_file else None
            drift_rows = [_predict_drift(str(p), settings, progress, args.only_step) for p in drift_paths]
        else:
            import concurrent.futures as cf
            with cf.ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
                drift_rows = list(pool.map(_predict_drift_isolated, [str(p) for p in drift_paths],
                                           [settings] * len(drift_paths)))
    drift_lines = []
    for r in drift_rows:
        if r.get("not_read") or r.get("failed"):
            what = f"not read: {r['not_read']}" if r.get("not_read") else r["failed"]
            drift_lines.append(f"| {r['spec']} | {r['field']} | {what} | | | |")
            continue
        t = r["targets"]
        drift_lines.append(f"| {r['spec']} | {r['field']} | {_drift_cell(t.get('symmetric-write', {}))} | "
                           f"{_drift_cell(t.get('threshold-write', {}))} | {_drift_cell(t.get('phase-locking', {}))} | "
                           f"{_completion_cell(r.get('completion', {}))} |")
    md = "# Predictions of the rules for the specifications of a directory\n\n"
    if lines:
        md += ("## Hysteron networks, the order of writes and the closure of quantum observables\n\n"
               "Quantum realizations (schema fieldbridge-quantum/1) are read with the closure rule: the closure of the "
               "observable, the index (Hamiltonian terms removed until it rotates), whether every line of the signal lies in "
               "the closure spectrum, and the spread of the lines (`quantum closure`). "
               "F counts the unsatisfiable edges of the network with the drive as a node. F = 0: the return to every turning "
               "point is exact. F > 0: states beyond return-point memory are possible (the rule is sufficient, not necessary; "
               "the one-dimensional antiferromagnetic chain of Deutsch, Dhar and Narayan keeps return-point memory with F > 0). "
               "The typical gain is a_k (F - F0_k) of the sequence-memory ensemble, a comparison value for random networks "
               "of similar thresholds and couplings only. Measured columns need the enumeration, feasible up to "
               f"{args.enumerate_up_to} elements here.\n\n"
               "| Specification | Elements | Mean degree | F | Balanced | Typical gain (bits) | Measured |\n"
               "| --- | --- | --- | --- | --- | --- | --- |\n" + "\n".join(lines) + "\n\n")
    if drift_lines:
        md += ("## Drift equations and networks: the write targets\n\n"
               "Each target is derived from the specification as `memory codiscover` derives it: reached (the model's write "
               "point carries the mechanism, with the law constant where the law was checked), or obstructed, with what "
               "obstructs it. Where the symmetric write is obstructed, the completion rung (`memory complete`) keeps what the "
               "obstruction names and derives the law of the completed model; its certificate is the full simulation. "
               f"Time limits: {args.budget} s per write target, {args.completion_budget} s for the completion"
               f"{' (not applied: the specifications ran in this process)' if args.no_isolate else ''}"
               f"{'; the law constants were not checked' if args.no_law else ''}.\n\n"
               "| Specification | Field | Symmetric write | Threshold write | Phase locking | Completion at the obstruction |\n"
               "| --- | --- | --- | --- | --- | --- |\n" + "\n".join(drift_lines) + "\n\n")
    if skipped:
        md += f"{skipped} specification(s) of other kinds were skipped (fields have their own command).\n\n"
    md += f"{BOUNDARY}\n"
    all_rows = rows + drift_rows
    _write(Path(args.out_dir), "predictions",
           {"command": "predictions", "rule": "R8", "rules": ["R8", "holonomy", "closure", "write point"], **_provenance(),
            "settings": {**settings, "isolated": not args.no_isolate, "workers": args.workers, "exact_up_to": args.exact_up_to, "enumerate_up_to": args.enumerate_up_to}, "rows": all_rows}, md)
    reached = {t: sum(str(r.get("targets", {}).get(t, {}).get("status", "")).startswith("reached") for r in drift_rows) for t in WRITE_TARGETS}
    completed = sum(1 for r in drift_rows for e in r.get("completion", {}).get("write_points", [])
                    if (e.get("law") or {}).get("certified") or any((u.get("law") or {}).get("certified") for u in e.get("unfoldings", [])))
    print(json.dumps({"out_dir": args.out_dir, "specifications": len(all_rows),
                      "balanced": sum(bool(r.get("balanced")) for r in rows), "commuting": sum(bool(r.get("commuting")) for r in rows),
                      "rotating": sum(bool(r.get("rotating")) for r in rows),
                      "measured": sum("gain" in r for r in rows), "drift": len(drift_rows), "reached": reached,
                      "completions_certified": completed, "skipped": skipped}, indent=2))
    return 0


def cmd_gallery(args) -> int:
    import numpy as np
    from . import compose, discovery, spec
    from . import hysterons as hy
    from .visual import hysteron_figure
    # every material under the examples folder, in subfolders too: drift equations and networks get a memory card,
    # hysterons a card of their return to a turning point; fields have their own command (memory field)
    reals, hysteron_specs, skipped = [], [], []
    root = Path(args.examples)
    for path in sorted(root.rglob("*.json")):
        rel = path.relative_to(root).as_posix()
        try:
            kind = json.loads(path.read_text(encoding="utf-8")).get("kind")
        except (json.JSONDecodeError, UnicodeDecodeError) as err:
            skipped.append(f"{rel}: {err}")
            continue
        if kind == "field":
            continue
        if kind == hy.KIND:
            hysteron_specs.append(path)
            continue
        try:
            reals.append(spec.load(path))
        except spec.SpecError as err:
            skipped.append(f"{rel}: {err}")
    loops = None if args.no_loops else compose.run_loops(np.random.default_rng(args.seed))
    out = Path(args.out_dir)
    cards = discovery.run(reals, out, seed=args.seed, figures=True, quick=args.quick, loops=loops)
    rows = [{"id": c["id"], "name": c["name"], "verdict": c["verdict"]} for c in cards]
    for k, path in enumerate(hysteron_specs):
        model = hy.from_spec(json.loads(path.read_text(encoding="utf-8")))
        c = hy.card(model, np.random.default_rng(args.seed), count=args.subloops)
        c["id"] = f"card{len(cards) + k:02d}"
        c["verdict"] = {"stores": f"{model.n} elements with two states each", "writes": "turning points of a slow drive",
                        "constructs": "return-point memory" if c["class"] == "return-point" else "return not exact",
                        "holds": "no thermal loss on the time of the drive (Law 3 neglected)",
                        "lock": c["check"]["verdict"]}
        (out / f"{c['id']}.json").write_text(json.dumps(c, indent=1, default=float), encoding="utf-8")
        hysteron_figure(c, out / f"{c['id']}.png")
        rows.append({"id": c["id"], "name": c["name"], "verdict": c["verdict"]})
        print(f"[{c['id']}] {model.name}: {c['verdict']['constructs']}", flush=True)
    report = {"command": "gallery", **_provenance(), "input": {"examples": str(args.examples), "quick": args.quick},
              "cards": rows, "skipped": skipped,
              "loops": None if loops is None else {"consistent": sum(r["consistent"] for r in loops), "total": len(loops)}}
    md = ("# Gallery of realizations\n\nOne card per example specification; `index.html` shows the memory cards in a "
          "browser, and the cards of hysterons are the figures of `memory hysterons`:\n\n"
          + "\n".join(f"- `{r['id']}` {r['name']}: {r['verdict']['constructs']}" for r in rows)
          + f"\n\n{BOUNDARY}\n")
    _write(out, "gallery", report, md)
    print(json.dumps({"out_dir": args.out_dir, "cards": len(rows), "skipped": len(skipped),
                      "page": str(out / "index.html")}, indent=2))
    return 0


def cmd_codiscover(args) -> int:
    from . import codiscovery, spec
    dirs = [args.examples] if args.examples else codiscovery.TARGET_INFO[args.target].get("examples",
                                                                                         ("examples/memory",))
    paths = [Path(p) for p in args.specs] or [p for d in dirs for p in sorted(Path(d).glob("*.json"))]
    reals, skipped = [], []
    for path in paths:
        try:
            reals.append(spec.load(path))
        except spec.SpecError as err:
            skipped.append(f"{path.name}: {err}")
    rep = codiscovery.codiscover(reals, args.seed, target=args.target,
                                 check_law=not args.no_law, n_traj=args.trajectories, replicates=args.replicates)
    out = Path(args.out_dir)
    report = {"command": "codiscover", **_provenance(),
              "input": {"specifications": [str(p) for p in paths], "target": args.target, "law": not args.no_law,
                        "trajectories": args.trajectories, "replicates": args.replicates, "seed": args.seed},
              "specifications_sha256": {r.name: r.spec_sha256 for r in reals}, **rep, "skipped": skipped}
    _write(out, "codiscover", report, codiscovery.markdown(rep) + f"\n{BOUNDARY}\n")
    if not args.no_figure:
        from .visual import codiscovery_figure
        codiscovery_figure(rep, out / "codiscover.png")
    s = rep["summary"]
    print(json.dumps({"out_dir": args.out_dir, "reached": s["reached"], "obstructed": s["obstructed"],
                      "fields_reached": s["fields_reached"], "derivation_classes": s["derivation_classes"],
                      **({"law_constant": [s["law_constant_mean"], s["law_constant_stderr"]]}
                         if "law_constant_mean" in s else {}), "skipped": len(skipped)}, indent=2))
    return 0


def cmd_new(args) -> int:
    from . import contribute
    out = Path(args.out or f"examples/memory/{args.name}.json")
    if out.exists() and not args.force:
        raise SystemExit(f"{out} exists; choose another name or pass --force")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(contribute.template(args.name, args.carrier), indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {out}. Replace the placeholders and the placeholder dynamics with your material's equations, then "
          f"run: python3 -B -m fieldbridge memory check {out}")
    return 0


def cmd_check(args) -> int:
    from . import contribute
    ok = True
    for path in args.specs:
        result = contribute.check(path, structure_only=args.structure_only, seed=args.seed)
        ok = ok and result["passed"]
        md = contribute.markdown(result, path)
        print(md)
        if args.out_dir:
            name = Path(path).stem
            report = {"command": "check", **_provenance(), "input": path, **result}
            _write(Path(args.out_dir), f"check_{name}", report, md)
    return 0 if ok else 1


def cmd_catalog(args) -> int:
    from . import contribute
    result = contribute.catalog(args.paths, seed=args.seed)
    md = contribute.catalog_markdown(result)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(md, encoding="utf-8")
    rows = sorted(result["materials"], key=lambda r: r["name"].lower())
    catalog_json = json.dumps({"materials": rows, "skipped": result["skipped"]}, indent=1) + "\n"
    out.with_suffix(".json").write_text(catalog_json, encoding="utf-8")
    print(json.dumps({"out": str(out), "materials": len(result["materials"]), "skipped": len(result["skipped"])},
                     indent=2))
    return 0


def add_parser(sub) -> None:
    memory = sub.add_parser("memory", help="Build, predict and transfer memory on any carrier (needs the memory extra).")
    msub = memory.add_subparsers(required=True)

    def common(p, spec_arg=True):
        if spec_arg:
            p.add_argument("spec", help="fieldbridge-memory/1 specification (JSON).")
        p.add_argument("--out-dir", required=True, help="Directory for the report and its provenance.")
        p.add_argument("--seed", type=int, default=20260923)

    p = msub.add_parser("predict", help="Predictions from structure alone; no simulation.")
    common(p)
    p.set_defaults(func=cmd_predict)
    p = msub.add_parser("card", help="Calculated memory card, with the structural predictions compared.")
    common(p)
    p.add_argument("--quick", action="store_true", help="Fewer trajectories; no swept-write check.")
    p.add_argument("--no-figure", action="store_true")
    p.set_defaults(func=cmd_card)
    p = msub.add_parser("attach", help="Detach a memory signature and attach it to another carrier.")
    p.add_argument("--from", dest="source", required=True, help="Source specification.")
    p.add_argument("--to", dest="target", required=True, help="Target specification.")
    common(p, spec_arg=False)
    p.set_defaults(func=cmd_attach)
    p = msub.add_parser("design", help="Solve for the setting of a material parameter that restores a symmetric write.")
    common(p)
    p.add_argument("--free", required=True, help="Material parameter to tune.")
    p.set_defaults(func=cmd_design)
    p = msub.add_parser("audit", help="Audit the symmetric-write law: each specific way it could be wrong (window of the "
                                      "fit, one sign of the bias, the element's own noise, time step, few quanta, "
                                      "published sign) settled by its own test.")
    common(p)
    p.add_argument("--trajectories", type=int, default=4000, help="Trajectories of each law measurement.")
    p.add_argument("--replicates", type=int, default=4, help="Independent runs into which they are split.")
    p.add_argument("--checks", nargs="+", default=None, help="Only these checks (by name).")
    p.set_defaults(func=cmd_audit)
    p = msub.add_parser("complete", help="The law at an obstruction: keep what the obstruction names and derive the "
                                         "write law of the completed model from the sweep, with its certificate.")
    common(p)
    p.add_argument("--rate", type=float, default=0.02, help="Rate at which the sweep drives the critical rate.")
    p.add_argument("--gamma", type=float, default=0.05,
                   help="Noise against the growth of the barrier, Gamma = D_s |a3| / rate (when --noise is not given).")
    p.add_argument("--noise", type=float, default=None, help="Noise intensity D (overrides --gamma).")
    p.add_argument("--floors", type=float, nargs="+", default=[0.5, 1.0, 2.0],
                   help="Bias strengths in units of the noise floor sqrt(D_s) rate^(1/4) / pi^(1/4).")
    p.add_argument("--samples", type=int, default=4000, help="Samples of the linear-stage Gaussian.")
    p.add_argument("--trajectories", type=int, default=2000, help="Trajectories of the certificate per bias.")
    p.add_argument("--no-certificate", action="store_true", help="Derive the law without the full simulation.")
    p.add_argument("--curvature", action="store_true", help="At a cusp, find the curvature of the sweep that cancels "
                   "its intrinsic offset (five curved sweeps and a check, each with the full dynamics).")
    p.set_defaults(func=cmd_complete)
    p = msub.add_parser("changes", help="Changes of a model (remove a term, make another parameter the control, extend "
                                        "the range of the control, add a term of a donor model) and pairs of changes, each judged by the symmetric "
                                        "and threshold writes; open pairs do what neither change does alone.")
    common(p)
    p.add_argument("--pairs", action="store_true", help="Also judge pairs of changes and their closure.")
    p.add_argument("--first", type=int, default=8, help="First changes of the pairs (steps first).")
    p.add_argument("--only", nargs="+", default=None,
                   help="Restrict the first changes of the pairs to those whose description contains one of these.")
    p.add_argument("--kinds", nargs="+", default=["remove", "control", "range"], choices=["remove", "control", "range"],
                   help="Kinds of change.")
    p.add_argument("--budget", type=int, default=60, help="Wall-clock seconds of each derivation in its child process.")
    p.add_argument("--donor", default=None, help="A second specification whose terms the changes of kind add take.")
    p.add_argument("--no-isolate", action="store_true", help="Derive in this process (no time limit; a crash in a "
                   "numerical routine ends the run).")
    p.add_argument("--save", action="store_true", help="Write the changed specifications of steps and open pairs.")
    p.add_argument("--verbose", action="store_true", help="Report each change as it is judged.")
    p.set_defaults(func=cmd_changes)
    p = msub.add_parser("phase", help="Memory in the phase of an oscillating realization: phase response, retention, writing.")
    common(p)
    p.add_argument("--kick", type=float, default=2.0, help="Kick on the first coordinate for the response curve.")
    p.add_argument("--shift", type=float, default=1.5707963, help="Written phase shift (rad).")
    p.add_argument("--noise", type=float, default=0.005, help="Noise intensity for the retention test.")
    p.set_defaults(func=cmd_phase)
    p = msub.add_parser("regimes", help="Information about a write for kappa > 0, = 0 and < 0 (exact and simulated).")
    common(p, spec_arg=False)
    p.add_argument("--kappa", type=float, default=1.0)
    p.add_argument("--s0", type=float, default=0.3, help="Distance of the written state from the reference.")
    p.add_argument("--noise", type=float, default=0.02)
    p.add_argument("--n", type=int, default=20000, help="Trajectories per regime.")
    p.set_defaults(func=cmd_regimes)
    p = msub.add_parser("noise-ratio", help="The swept-write law across the noise ratio Gamma = D|a3|/r: pi^(1/4) "
                                            "for weak noise, an activated law for strong noise; continuous and "
                                            "first-order switches.")
    p.add_argument("spec", nargs="?", default=None, help="fieldbridge-memory/1 specification (JSON).")
    p.add_argument("--out-dir", required=True, help="Directory for the report and its provenance.")
    p.add_argument("--seed", type=int, default=20260923)
    p.add_argument("--law", action="store_true", help="Print the law c(Gamma) without a specification.")
    p.add_argument("--gammas", default="0.05,1,10", help="Noise ratios, comma separated.")
    p.add_argument("--q", type=float, default=5e-4,
                   help="D_s|a3|, small enough that the choice is made near the bifurcation (default 5e-4).")
    p.add_argument("--trajectories", type=int, default=400, help="Sweeps per side of the bias.")
    p.set_defaults(func=cmd_noise_ratio)
    p = msub.add_parser("dynamics", help="Chaotic loss of the initial state (predictability law) and blow-up in "
                                         "finite time (exponent from the dominant balance).")
    common(p)
    p.add_argument("--starts", type=int, default=6, help="Random starts for each test.")
    p.set_defaults(func=cmd_dynamics)
    p = msub.add_parser("field", help="Memory in a field: conservation and dimension set the law of loss.")
    common(p, spec_arg=False)
    p.add_argument("spec", help="Field specification (examples/memory/fields/*.json).")
    p.add_argument("--pairs", type=int, default=120, help="Written/unwritten pairs in the nonlinear simulation.")
    p.set_defaults(func=cmd_field)
    p = msub.add_parser("hysterons", help="Return to a turning point: interacting hysterons under a slow drive, "
                                          "the structural prediction against subloops.")
    common(p, spec_arg=False)
    p.add_argument("spec", help="Specification of kind hysterons (examples/memory/hysterons/*.json), or of kind "
                                "equations with a block slow_drive (examples/slow_drive/*.json).")
    p.add_argument("--subloops", type=int, default=24, help="Subloops for each order of relaxation, simple and nested.")
    p.add_argument("--no-figure", action="store_true", help="Skip the figure hysterons.png.")
    p.set_defaults(func=cmd_hysterons)
    p = msub.add_parser("counterpart", help="The balanced counterpart of a hysteron network: the smallest reversal of "
                                            "signs that guarantees the return to every turning point, and the "
                                            "states of the network measured against it.")
    common(p, spec_arg=False)
    p.add_argument("spec", help="Specification of kind hysterons (examples/counterpart/*.json).")
    p.add_argument("--max-states", type=int, default=20000,
                   help="States of a transition graph before its enumeration stops (the result then says so).")
    p.add_argument("--exact-up-to", type=int, default=16,
                   help="Elements up to which the frustration index is exact; above, an upper bound by local search.")
    p.add_argument("--no-return-test", action="store_true", help="Skip the exhaustive return test.")
    p.set_defaults(func=cmd_counterpart)

    p = msub.add_parser("order", help="The order of writes in a bistable network measured against its uncoupled "
                                      "reference: the holonomy rule as a measurement.")
    common(p, spec_arg=False)
    p.add_argument("spec", help="Specification of kind order (examples/counterpart/order/*.json).")
    p.add_argument("--exact-up-to", type=int, default=10, help="Nodes up to which the commutation index is exact.")
    p.set_defaults(func=cmd_order)
    p = msub.add_parser("record", help="The record of every drift specification of a file or a directory: the write "
                                       "targets, the write points with their normal forms, the completions with their "
                                       "classes and retention, and the motion of a realization without a stable state "
                                       "(schema fieldbridge-record/1); every step in a child process under a time limit.")
    p.add_argument("path", nargs="?", default=None, help="A specification file or a directory of specifications "
                   "(subfolders with --recursive); not needed with --from-records.")
    p.add_argument("--out-dir", default="out/memory_record")
    p.add_argument("--recursive", action="store_true", help="Include the specifications of the subfolders.")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--budget", type=int, default=300, help="Wall-clock seconds per write target and for the write points.")
    p.add_argument("--completion-budget", type=int, default=1200, help="Wall-clock seconds for the completion rung.")
    p.add_argument("--attractor-budget", type=int, default=900, help="Wall-clock seconds for the motion of a realization "
                   "without a stable state.")
    p.add_argument("--no-law", action="store_true", help="Derive the write targets without checking their law constants.")
    p.add_argument("--trajectories", type=int, default=400, help="Trajectories of the law check per target.")
    p.add_argument("--samples", type=int, default=4000, help="Samples of the completion's linear-stage Gaussian.")
    p.add_argument("--trajectories-completion", type=int, default=2000, help="Trajectories of the completion's certificate per bias.")
    p.add_argument("--workers", type=int, default=1, help="Specifications recorded at a time (each step still in its own process).")
    p.add_argument("--no-isolate", action="store_true", help="Run the steps in this process (no time limits).")
    p.add_argument("--steps", default=None, help="Comma-separated steps to run (default all: "
                   + ",".join(STEPS_FOR_PARSER) + "); the others are recorded as not run.")
    p.add_argument("--from-records", default=None, metavar="FILE", help="No derivation: write the table of the write "
                   "targets and completions (<stem>.md, <stem>.json in --out-dir) from an existing records.jsonl or records.json.")
    p.add_argument("--intro", default=None, help="With --from-records: a sentence placed before the counts of the table.")
    p.add_argument("--only-step", default=None, choices=list(STEPS_FOR_PARSER), help=argparse.SUPPRESS)
    p.set_defaults(func=cmd_record)
    p = msub.add_parser("predictions", help="What the rules predict for every specification in a directory: for "
                                           "hysterons rule R8 (frustration index with the drive, balance, the typical "
                                           "gain of the sequence-memory ensemble, the measured gain and return test "
                                           "where the enumeration is feasible); for bistable networks the holonomy of "
                                           "the writes; for drift equations and networks the three write targets under "
                                           "a time budget and the completion rung at an obstructed symmetric write.")
    p.add_argument("directory", help="Directory of specifications (e.g. examples/memory/hysterons, examples/counterpart, "
                                      "examples/quantum/semiclassical); its subfolders only with --recursive.")
    p.add_argument("--recursive", action="store_true", help="Include the specifications of the subfolders.")
    p.add_argument("--out-dir", default="out/memory_predictions")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--exact-up-to", type=int, default=16, help="Elements up to which the frustration index is exact.")
    p.add_argument("--enumerate-up-to", type=int, default=16, help="Elements up to which the states are enumerated.")
    p.add_argument("--max-states", type=int, default=20000)
    p.add_argument("--budget", type=int, default=300, help="Wall-clock seconds allowed to the child process that "
                   "derives one write target of one drift specification (0: no limit).")
    p.add_argument("--completion-budget", type=int, default=900, help="Wall-clock seconds allowed to the child "
                   "process that runs the completion rung on one specification.")
    p.add_argument("--no-law", action="store_true", help="Derive the write targets without checking their law constants.")
    p.add_argument("--trajectories", type=int, default=400, help="Trajectories of the law check per target.")
    p.add_argument("--no-complete", action="store_true", help="Skip the completion rung at an obstructed symmetric write.")
    p.add_argument("--samples", type=int, default=4000, help="Samples of the completion's linear-stage Gaussian.")
    p.add_argument("--trajectories-completion", type=int, default=2000, help="Trajectories of the completion's certificate per bias.")
    p.add_argument("--workers", type=int, default=1, help="Child processes over the drift specifications.")
    p.add_argument("--no-isolate", action="store_true", help="Derive the drift specifications in this process instead "
                   "of one child process each (no time limit applies then, and a crash in a numerical routine ends "
                   "the run).")
    p.add_argument("--progress-file", default=None, help=argparse.SUPPRESS)
    p.add_argument("--only-step", default=None, choices=list(WRITE_TARGETS) + ["completion"], help=argparse.SUPPRESS)
    p.set_defaults(func=cmd_predictions)
    p = msub.add_parser("gallery", help="Evaluate every example specification and write an HTML gallery with figures.")
    common(p, spec_arg=False)
    p.add_argument("--examples", default="examples/memory",
                   help="Directory of specifications, read with its subfolders: equations and networks get a memory "
                        "card, hysterons the card of memory hysterons; fields are left to memory field.")
    p.add_argument("--quick", action="store_true", help="Fewer trajectories; no swept-write check.")
    p.add_argument("--subloops", type=int, default=24, help="Subloops of each card of hysterons, per procedure.")
    p.add_argument("--no-loops", action="store_true", help="Omit the loops section.")
    p.set_defaults(func=cmd_gallery)
    p = msub.add_parser("loops", help="Loops of genes, spins and rotors: holonomy against simulation.")
    common(p, spec_arg=False)
    p.add_argument("--no-figure", action="store_true")
    p.set_defaults(func=cmd_loops)
    p = msub.add_parser("networks", help="Random networks: structural predictions against simulation.")
    common(p, spec_arg=False)
    p.add_argument("--count", type=int, default=40, help="Networks per kind.")
    p.set_defaults(func=cmd_networks)
    p = msub.add_parser("codiscover", help="Derive one mechanism in realizations from different fields and compare "
                                           "the derivations and their invariants.")
    common(p, spec_arg=False)
    p.add_argument("specs", nargs="*", help="Specifications (default: every file in --examples).")
    p.add_argument("--examples", default=None,
                   help="Directory of specifications used when none is given (default: examples/memory; for "
                        "phase locking also examples/memory/oscillators).")
    p.add_argument("--target", default="symmetric-write", choices=["symmetric-write", "threshold-write", "phase-locking"],
                   help="Mechanism to derive: a symmetric write (pitchfork), a threshold write (fold) or phase "
                        "locking (Adler equation).")
    p.add_argument("--trajectories", type=int, default=400,
                   help="Trajectories per realization for the swept-write law (symmetric write only).")
    p.add_argument("--replicates", type=int, default=8,
                   help="Independent runs into which the trajectories are split; the error of each law constant is "
                        "the larger of their scatter and the binomial error (symmetric write only).")
    p.add_argument("--no-law", action="store_true", help="Derivations and canonical forms only; no simulation.")
    p.add_argument("--no-figure", action="store_true")
    p.set_defaults(func=cmd_codiscover)
    p = msub.add_parser("new", help="Write a template specification for a new material.")
    p.add_argument("name", help="Short name of the material; also the file name (examples/memory/NAME.json).")
    p.add_argument("--carrier", default="euclid", choices=["euclid", "orthant", "torus"],
                   help="Kind of state variables: real numbers, non-negative amounts, or angles.")
    p.add_argument("--out", help="Output file (default examples/memory/NAME.json).")
    p.add_argument("--force", action="store_true", help="Overwrite an existing file.")
    p.set_defaults(func=cmd_new)
    p = msub.add_parser("check", help="Check that a specification is ready to contribute; exit status 1 if not.")
    p.add_argument("specs", nargs="+", help="Specification files.")
    p.add_argument("--structure-only", action="store_true", help="Skip the quick memory card.")
    p.add_argument("--out-dir", help="Also write check_<name>.json and .md here.")
    p.add_argument("--seed", type=int, default=20260923)
    p.set_defaults(func=cmd_check)
    p = msub.add_parser("catalog", help="Write the table of materials (docs/materials.md and .json) from the "
                                        "specifications.")
    p.add_argument("paths", nargs="*", default=["examples/memory"], help="Specification files or directories.")
    p.add_argument("--out", default="docs/materials.md")
    p.add_argument("--seed", type=int, default=20260923)
    p.set_defaults(func=cmd_catalog)
