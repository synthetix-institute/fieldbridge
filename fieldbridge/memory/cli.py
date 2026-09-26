"""``fieldbridge memory ...``: build, predict and transfer memory from a specification file.

  predict  SPEC                structure only: what the realization can store and how (no simulation)
  card     SPEC                the calculated memory card, the structural predictions and their comparison
  attach   --from SPEC --to SPEC   which properties of a memory are preserved, changed or obstructed on another carrier
  design   SPEC --free PARAM   the setting of PARAM that restores a symmetric (pitchfork) write
  phase    SPEC                memory in the phase of an oscillating realization (phase response, retention, writing)
  regimes                      information about a write for kappa > 0, = 0 and < 0
  field    SPEC                memory in a field: conservation and dimension set the law of loss
  gallery                      every example: cards, figures and an HTML gallery (index.html)
  loops                        loops of genes, spins and rotors: holonomy against simulation
  networks                     random networks of genes, spins and rotors: structure against simulation
  codiscover [SPEC ...]        one mechanism (--target symmetric-write, threshold-write or phase-locking) derived in
                               realizations from different fields: derivations, obstructions, invariants (default:
                               every example)

Contributing a material (CONTRIBUTING.md):
  new NAME                     a template specification that loads and runs, with placeholders to replace
  check SPEC [SPEC ...]        required and recommended checks, the structural predictions and a quick card
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
        lines.append(f"- no setting of {free} removes the obstruction within the search range")
    else:
        d = c["direction"]
        cstate = ", ".join(f"{x:.4g}" for x in c["state"])
        lines += [f"- the obstruction vanishes at {free} = {c['value']:.4g}, {control} = {c['control']:.4g}, "
                  f"state ({cstate})",
                  f"- sweep through this point along ({control}, {free}) proportional to "
                  f"({d[control]:.3g}, {d[free]:.3g})",
                  f"- kind along the sweep: {c['kind_along_sweep']}"]
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


def cmd_gallery(args) -> int:
    import numpy as np
    from . import compose, discovery, spec
    reals, skipped = [], []
    for path in sorted(Path(args.examples).glob("*.json")):
        try:
            reals.append(spec.load(path))
        except spec.SpecError as err:
            skipped.append(f"{path.name}: {err}")
    loops = None if args.no_loops else compose.run_loops(np.random.default_rng(args.seed))
    out = Path(args.out_dir)
    cards = discovery.run(reals, out, seed=args.seed, figures=True, quick=args.quick, loops=loops)
    rows = [{"id": c["id"], "name": c["name"], "verdict": c["verdict"]} for c in cards]
    report = {"command": "gallery", **_provenance(), "input": {"examples": str(args.examples), "quick": args.quick},
              "cards": rows, "skipped": skipped,
              "loops": None if loops is None else {"consistent": sum(r["consistent"] for r in loops), "total": len(loops)}}
    md = ("# Gallery of realizations\n\nOpen `index.html` in a browser. One card per example specification:\n\n"
          + "\n".join(f"- `{r['id']}` {r['name']}: {r['verdict']['constructs']}" for r in rows)
          + f"\n\n{BOUNDARY}\n")
    _write(out, "gallery", report, md)
    print(json.dumps({"out_dir": args.out_dir, "cards": len(rows), "skipped": len(skipped),
                      "page": str(out / "index.html")}, indent=2))
    return 0


def cmd_codiscover(args) -> int:
    import numpy as np
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
    rep = codiscovery.codiscover(reals, np.random.default_rng(args.seed), target=args.target,
                                 check_law=not args.no_law, n_traj=args.trajectories)
    out = Path(args.out_dir)
    report = {"command": "codiscover", **_provenance(),
              "input": {"specifications": [str(p) for p in paths], "target": args.target, "law": not args.no_law,
                        "trajectories": args.trajectories, "seed": args.seed},
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
    p = msub.add_parser("field", help="Memory in a field: conservation and dimension set the law of loss.")
    common(p, spec_arg=False)
    p.add_argument("spec", help="Field specification (examples/memory/fields/*.json).")
    p.add_argument("--pairs", type=int, default=120, help="Written/unwritten pairs in the nonlinear simulation.")
    p.set_defaults(func=cmd_field)
    p = msub.add_parser("gallery", help="Evaluate every example specification and write an HTML gallery with figures.")
    common(p, spec_arg=False)
    p.add_argument("--examples", default="examples/memory", help="Directory of specifications (kinds equations, network).")
    p.add_argument("--quick", action="store_true", help="Fewer trajectories; no swept-write check.")
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
