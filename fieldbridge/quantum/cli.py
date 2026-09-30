"""``fieldbridge quantum ...``: the language of mechanisms on quantum carriers.

  detach SPEC                        the carrier-free signature of a realization and what stays with its carrier
  attach --from SPEC --to CARRIER    write the detached mechanism on another carrier, in that carrier's operators,
                                     and derive it there again
  carriers                           the carriers available to attach
  codiscover [SPEC ...]              derive the Bloch rotation in realizations from different fields: derivations
                                     through the algebra or through the closure of the observable, obstructions and
                                     invariants (default: every file in examples/quantum)

Every command writes <name>.json and <name>.md to --out-dir, with the input hashes and the hash of the implementation.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Dict, List

BOUNDARY = ("Exact calculations for the supplied Hamiltonians on finite carriers. They do not establish that a model "
            "describes a physical system, and they do not establish novelty.")


def _provenance(specs: List[Dict] = ()) -> Dict[str, object]:
    here = Path(__file__).resolve().parent
    impl = hashlib.sha256(b"".join(p.read_bytes() for p in sorted(here.glob("*.py")))).hexdigest()
    versions = {"python": platform.python_version()}
    for mod in ("numpy", "sympy", "matplotlib"):
        try:
            versions[mod] = __import__(mod).__version__
        except ImportError:
            versions[mod] = None
    hashes = {s["name"]: hashlib.sha256(json.dumps(s, sort_keys=True).encode()).hexdigest() for s in specs}
    return {"implementation_sha256": impl, "versions": versions, "input_sha256": hashes, "novelty_established": False,
            "evidence_boundary": BOUNDARY}


def _write(out: Path, name: str, report: Dict, md: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.json").write_text(json.dumps(report, indent=2, default=float) + "\n", encoding="utf-8")
    (out / f"{name}.md").write_text(md, encoding="utf-8")


def _terms_text(terms: List[List]) -> str:
    return " ".join(f"{'+' if c >= 0 else '-'} {abs(c):.4g} {name}" for c, name in terms).lstrip("+ ")


def _strip(row: Dict) -> Dict:
    """A derivation without the time series of the law (kept in the figure, not in the report)."""
    law = {k: v for k, v in row.get("law", {}).items() if k not in ("times", "f_exact", "f_law")}
    return {**row, "law": law}


def cmd_detach(args) -> int:
    from . import language as ql
    real = ql.load(args.spec)
    d = ql.detach(real)
    report = {"command": "quantum detach", **_provenance([real.spec]),
              **{k: v for k, v in d.items() if k != "derivation"}, "derivation": _strip(d["derivation"])}
    if not d["detached"]:
        md = f"# Detachment: {real.name}\n\nNothing to detach: {d['reason']}.\n\n{BOUNDARY}\n"
    else:
        s, lb, row = d["signature"], d["left_behind"], d["derivation"]
        if row["status"] == ql.REACHED:  # through the algebra: H and the observable generate su(2)
            how = ""
            relations, angle = f"algebra su(2): {s['relations']}", "angle between Omega and the observable"
            left = f"- representation: {ql._rep(lb['representation'])}"
            if lb.get("closing_operators"):
                left += ("\n- operators that close under commutation (A): " + ", ".join(lb["closing_operators"])
                         + " (each divided by 2 is a spin-1/2 component)")
        else:
            how = "; the rotation is reached through the closure of the observable"
            relations = f"rotation of the operators of the closure: {s['relations']}"
            angle = "angle between the rotation axis and the observable"
            left = f"- {row['algebra']}"
            if lb.get("algebra_operators"):
                left += "\n- operators of this algebra (A): " + ", ".join(lb["algebra_operators"])
            if lb.get("closure_operators"):
                left += "\n- closure of the observable (O): " + ", ".join(lb["closure_operators"])
            if lb["representation"]:
                left += f"\n- the closure is an su(2); representation: {ql._rep(lb['representation'])}"
        if "generators" in lb:
            left += ("\n- the rotating vector in the canonical frame (K; J3 along Omega):\n"
                     + "\n".join(f"  - J{a + 1} = {_terms_text(t)}" for a, t in enumerate(lb["generators"])))
        md = (f"# Detachment: {real.name}\n\n{real.spec['question']}\n\n"
              f"Derivation: {row['word']}{how}.\n\n"
              f"## Signature (independent of the carrier)\n\n"
              f"- {relations}\n- rotation rate |Omega| = {s['rate']:.6g}\n"
              f"- weight of the observable |n| = {s['weight']:.6g}\n"
              f"- {angle}: {s['theta_deg']:.4g} degrees\n"
              f"- law: {s['law']}; the observable is inverted as far as it can be at t = pi/|Omega| = "
              f"{d['law']['inversion_time']:.6g}; exact evolution deviates by {d['law']['residual']:.1e}\n\n"
              f"## Left with the carrier\n\n- carrier: {lb['carrier']} ({lb['field']})\n"
              f"- Hilbert space dimension {lb['hilbert_dimension']}; sector dimension {lb['sector_dimension']}\n"
              f"{left}\n\n{BOUNDARY}\n")
    _write(Path(args.out_dir), "detach", report, md)
    print(json.dumps({"out_dir": args.out_dir, "detached": d["detached"],
                      **({"signature": d["signature"]} if d["detached"] else {"reason": d["reason"]})}, indent=2))
    return 0


def cmd_attach(args) -> int:
    from . import language as ql
    src = ql.load(args.source)
    a = ql.attach(src, args.to, args.size)
    out = Path(args.out_dir)
    if not a["attached"] and "spec" not in a:
        report = {"command": "quantum attach", **_provenance([src.spec]), **a}
        _write(out, "attach", report, f"# Attachment\n\nNothing to attach: {a['reason']}.\n\n{BOUNDARY}\n")
        print(json.dumps({"out_dir": args.out_dir, "attached": False, "reason": a["reason"]}, indent=2))
        return 1
    out.mkdir(parents=True, exist_ok=True)
    (out / "attached.json").write_text(json.dumps(a["spec"], indent=2) + "\n", encoding="utf-8")
    report = {"command": "quantum attach", **_provenance([src.spec, a["spec"]]),
              **{k: v for k, v in a.items() if k not in ("source_derivation", "target_derivation")},
              "source_derivation": _strip(a["source_derivation"]), "target_derivation": _strip(a["target_derivation"])}
    p, c = a["preserved"], a["changed"]
    native = "\n".join(f"- {k.replace('_', ' ')}: {v}" for k, v in a["native"].items())
    through = ""
    if a["source_derivation"]["status"] == ql.REACHED_CLOSURE:
        through = (f" In the source it is reached through the closure of the observable "
                   f"({a['source_derivation']['word']}); on the new carrier the Hamiltonian and the observable generate "
                   f"su(2).")
    md = (f"# Attachment: {a['source']} -> {ql.CARRIERS[a['to']]}\n\n"
          f"The rotation detached from '{a['source']}' is written in the operators of the new carrier and derived there "
          f"again ({a['target_derivation']['word']}).{through}\n\n## Written in the carrier's own operators\n\n{native}\n\n"
          f"The specification is saved as attached.json.\n\n## Preserved and changed\n\n"
          f"| property | source | target |\n|---|---|---|\n"
          f"| rotation rate | {p['rate']['source']:.6g} | {p['rate']['target']:.6g} |\n"
          f"| weight of the observable | {p['weight']['source']:.6g} | {p['weight']['target']:.6g} |\n"
          f"| angle (degrees) | {p['theta_deg']['source']:.6g} | {p['theta_deg']['target']:.6g} |\n"
          f"| law, deviation of exact evolution | {a['law']['source_residual']:.1e} | {a['law']['target_residual']:.1e} |\n"
          f"| carrier | {c['carrier']['source']} | {c['carrier']['target']} |\n"
          f"| dimension of the space it acts on | {c['sector_dimension']['source']} | {c['sector_dimension']['target']} |\n"
          f"| representation | {ql._rep(c['representation']['source'])} | {ql._rep(c['representation']['target'])} |\n\n"
          f"{BOUNDARY}\n")
    _write(out, "attach", report, md)
    if not args.no_figure:
        try:
            from .figure import attach_figure
            attach_figure(a, out / "attach.png")
        except ImportError:
            pass
    print(json.dumps({"out_dir": args.out_dir, "attached": a["attached"], "native": a["native"],
                      "law_residual": a["law"]["target_residual"]}, indent=2))
    return 0


def cmd_carriers(args) -> int:
    from .language import CARRIERS
    print("\n".join(f"{k:16s} {v}" for k, v in CARRIERS.items()))
    return 0


def cmd_codiscover(args) -> int:
    from . import language as ql
    paths = [Path(p) for p in args.specs] or sorted(Path(args.examples).glob("*.json"))
    reals, skipped = [], []
    for p in paths:
        try:
            reals.append(ql.load(p))
        except ql.SpecError as err:
            skipped.append(f"{p.name}: {err}")
    rep = ql.codiscover(reals)
    out = Path(args.out_dir)
    report = {"command": "quantum codiscover", **_provenance([r.spec for r in reals]),
              "summary": rep["summary"], "rows": [_strip(r) for r in rep["rows"]], "skipped": skipped}
    _write(out, "codiscover", report, ql.markdown(rep) + f"\n{BOUNDARY}\n")
    if not args.no_figure:
        try:
            from .figure import codiscovery_figure
            codiscovery_figure(rep, out / "codiscover.png")
        except ImportError:
            pass
    s = rep["summary"]
    print(json.dumps({"out_dir": args.out_dir, "reached": s["reached"],
                      "reached_through_closure": s["reached_through_closure"], "obstructed": s["obstructed"],
                      "fields_reached": s["fields_reached"], "law_residual_max": s["law_residual_max"],
                      "skipped": len(skipped)}, indent=2))
    return 0


def add_parser(sub) -> None:
    import numpy  # noqa: F401  (the quantum language needs numpy and sympy: pip install -e '.[construction]')
    import sympy  # noqa: F401
    q = sub.add_parser("quantum", help="The language of mechanisms on quantum carriers: detach, attach, co-discover.")
    qsub = q.add_subparsers(required=True)
    p = qsub.add_parser("detach", help="The carrier-free signature of a realization and what stays with its carrier.")
    p.add_argument("spec", help="fieldbridge-quantum/1 specification (JSON).")
    p.add_argument("--out-dir", required=True)
    p.set_defaults(func=cmd_detach)
    p = qsub.add_parser("attach", help="Write a detached mechanism on another carrier and derive it there.")
    p.add_argument("--from", dest="source", required=True, help="Source specification.")
    p.add_argument("--to", required=True, help="Carrier: see 'fieldbridge quantum carriers'.")
    p.add_argument("--size", type=int, help="Number of spins, particles or 2j, where the carrier needs one.")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--no-figure", action="store_true")
    p.set_defaults(func=cmd_attach)
    p = qsub.add_parser("carriers", help="List the carriers available to attach.")
    p.set_defaults(func=cmd_carriers)
    p = qsub.add_parser("codiscover", help="Derive the Bloch rotation in realizations from different fields.")
    p.add_argument("specs", nargs="*", help="Specifications (default: every file in --examples).")
    p.add_argument("--examples", default="examples/quantum")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--no-figure", action="store_true")
    p.set_defaults(func=cmd_codiscover)
