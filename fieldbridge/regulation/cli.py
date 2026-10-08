"""``fieldbridge regulation``: cards, surveys and checks of certificates for regulated realizations.

    fieldbridge regulation card examples/regulation/chemotaxis_tu2008.json --out-dir out/chemotaxis
    fieldbridge regulation survey examples/regulation --out-dir out/regulation
    fieldbridge regulation check out/chemotaxis/regulation.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

from ..core.provenance import provenance, write_report

BOUNDARY = ("Calculated and structural statements about a supplied model. They do not establish that the model "
            "describes a physical system, and they do not establish novelty; both need comparison with experiment "
            "and with the literature.")
HERE = Path(__file__).resolve().parent


def _plain(x):
    import numpy as np
    if isinstance(x, dict):
        return {k: _plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.ndarray):
        return _plain(x.tolist())
    if isinstance(x, complex):
        return [x.real, x.imag]
    return x


def _phi_text(integ: Dict) -> str:
    c = integ["coefficients"]
    terms = [f"{v:+.4g} {k}" for k, v in c["w"].items()] + [f"{v:+.4g} ln {k}" for k, v in c["v"].items()]
    phi = " ".join(terms).lstrip("+")
    if integ["stage"] == 1:
        return f"phi = {phi}, dphi/dt = {integ['k_I']:.4g} (y - {integ['set_point']:.6g})"
    if integ["stage"] == 2:
        g = integ["gain_coefficients"]
        top = max(abs(v) for v in g.values())
        gain = " ".join(f"{v:+.4g}{'' if k == '1' else ' ' + k}" for k, v in g.items() if abs(v) > 1e-9 * top)
        return f"phi = {phi}, dphi/dt = ({gain.lstrip('+')}) (y - {integ['set_point']:.6g})"
    return f"phi = {phi}, dphi/dt = psi(y) with psi({integ['set_point']:.6g}) = 0, of the sign of y - y0"


def markdown(res: Dict) -> str:
    lines = [f"# Return to a set point: {res['name']}", "", f"Field: {res['field']}. Source: {res['source']}", "",
             f"Input: `{res['input']}` stepped from {res['u0']:.6g} to {', '.join(f'{u:.6g}' for u in res['steps'][1:])}. "
             f"Output: `{res['output']}`.", "", f"## Class: {res['class']}", ""]
    if res.get("gain_ratio") is not None:
        lines += [f"Static gain of the output G = {res['G']:.6g}; against the open-loop gain of the clamps, "
                  f"G/G_open = {res['gain_ratio']:.4g}." + (f" Relative sensitivity S = (u/y) G = {res['S']:.4g}."
                                                            if res.get("S") is not None else ""), ""]
    integ = res.get("integrator") or {}
    if integ.get("found"):
        lines += ["## Integrator", "", f"Stage {integ['stage']} (gain {integ['gain']}): {_phi_text(integ)}.", ""]
    elif res["class"] in ("perfect-adaptation", "fine-tuned-adaptation", "partial-adaptation", "no-adaptation"):
        lines += ["## Integrator", "", "None in the searched coordinates (linear and logarithmic, a gain linear in the "
                  "state, a rate that depends on the state only through the output).", ""]
    if res.get("note"):
        lines += [f"Note: {res['note']}.", ""]
    att = [a for a in res.get("attenuation") or []
           if a.get("remaining_fraction") is not None or a.get("open_loop_zero")]
    if att:
        lines += ["## Clamps", "", "| Variable | Role | Remaining fraction G/G_open |", "| --- | --- | --- |"]
        for a in att:
            fraction = ("undefined: G_open = 0, no response with the variable clamped" if a.get("open_loop_zero")
                        else f"{a['remaining_fraction']:.4g}")
            lines.append(f"| {a['variable']} | {a['role']} | {fraction} |")
        lines.append("")
    rob = res.get("robustness")
    if rob:
        lines += [f"Parameter points (each free parameter multiplied by a factor in [1/2, 2]): {rob['stable']} of "
                  f"{rob['samples']} with a stable steady state; G = 0 at {rob['gain_zero']} of them; an integrator at "
                  f"{rob['integrator_found']}.", ""]
    resp = [r for r in res.get("step_responses") or [] if "peak" in r]
    if resp:
        lines += ["## Steps", "", "| Input after the step | Peak deviation | Final deviation | Final / peak | Return time "
                  "| Calibration (integral / change of phi) |", "| --- | --- | --- | --- | --- | --- |"]
        for r in resp:
            cal = r.get("calibration_ratio")
            fp = r.get("final_over_peak")
            lines.append(f"| {r['u1']:.6g} | {r['peak']:.4g} | {r['final']:.4g} | {'' if fp is None else f'{fp:.4g}'} | "
                         f"{r['return_time']:.4g} | {'' if cal is None else f'{cal:.12f}'} |")
        lines.append("")
    lines += ["## Boundary", "", BOUNDARY, ""]
    return "\n".join(lines)


def _report(res: Dict, real) -> Dict:
    from .certificate import certificate
    rep = {"command": "regulation card", **provenance(HERE, BOUNDARY, real.spec), "card": _plain(res),
           "certificate": _plain(certificate(res))}
    return rep


def cmd_card(args) -> int:
    from .card import card
    from .spec import load
    real = load(args.spec)
    res = card(real, samples=args.samples, seed=args.seed)
    write_report(Path(args.out_dir), "regulation", _report(res, real), markdown(res))
    print(f"{real.name}: {res['class']} -> {Path(args.out_dir) / 'regulation.md'}")
    return 0


def _num(x) -> str:
    """A number for a table: values below 1e-9 in magnitude are zero to the accuracy of the calculation."""
    return "" if x is None else ("0" if abs(x) < 1e-9 else f"{x:.3g}")


def _integer_multiple(values: List[float], largest: int = 12):
    """The smallest m <= largest for which every m*v is an integer to 1e-6, or None."""
    for mult in range(1, largest + 1):
        if all(abs(mult * v - round(mult * v)) < 1e-6 for v in values):
            return mult
    return None


def _short_phi(integ: Dict) -> str:
    if not integ or not integ.get("found"):
        return "none found"
    c = integ["coefficients"]
    coef = list(c["w"].values()) + list(c["v"].values())
    mult = _integer_multiple(coef) or 1
    fmt = (lambda v: f"{round(mult * v):+d}".replace("+1 ", "+").replace("-1 ", "-")) if _integer_multiple(coef) \
        else (lambda v: f"{v:+.3g}")
    terms = [f"{fmt(v)} {k}" for k, v in c["w"].items()] + [f"{fmt(v)} ln {k}" for k, v in c["v"].items()]
    text = " ".join(terms).replace("+1 ", "+").replace("-1 ", "-").lstrip("+")
    if integ["stage"] == 1:
        gain = f"constant gain {mult * integ['k_I']:.3g}" if mult != 1 else f"constant gain {integ['k_I']:.3g}"
    elif integ["stage"] == 2:
        g = {k: v * mult for k, v in integ["gain_coefficients"].items() if abs(v) > 1e-9}
        gain = "gain " + " ".join(f"{v:+.3g}{'' if k == '1' else ' ' + k}" for k, v in g.items()).lstrip("+")
    else:
        gain = "rate a function of the output"
    return f"{text} ({gain})"


def tables_from_json(path: Path) -> List[str]:
    """The tables of a survey written earlier (the specifications are read again for their regulation blocks)."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    rows = [{"report": r, "spec_regulation": json.loads(Path(r["card"]["spec"]).read_text(encoding="utf-8"))
             ["regulation"]} for r in data["results"]]
    return _survey_tables(rows)


def _survey_tables(rows: List[Dict]) -> List[str]:
    """Two tables: the published models, and the controls with the change each one makes."""
    models = ["| Model | Field | Class | Integrator | G/G_open | Final / peak after the first step |",
              "| --- | --- | --- | --- | --- | --- |"]
    controls = ["| Control | Change | Class |", "| --- | --- | --- |"]
    for r in rows:
        c = r["report"]["card"]
        reg = r["spec_regulation"]
        ratio = c.get("gain_ratio")
        resp = [s for s in c.get("step_responses") or [] if s.get("final_over_peak") is not None]
        fp = _num(resp[0]["final_over_peak"]) if resp else ""
        if reg.get("control_of"):
            controls.append(f"| {c['name']} | {reg.get('change', '')} | {c['class']} |")
        else:
            models.append(f"| {c['name']} | {c['field']} | {c['class']} | {_short_phi(c.get('integrator'))} | "
                          f"{_num(ratio)} | {fp} |")
    return ["# Regulation survey", "", "## Models", ""] + models + ["", "## Controls", ""] + controls


def cmd_survey(args) -> int:
    from .card import card, summary
    from .spec import load
    files: List[Path] = []
    for item in args.paths:
        path = Path(item)
        files += sorted(path.glob("*.json")) if path.is_dir() else [path]
    rows, bad = [], []
    for f in files:
        real = load(f)
        res = card(real, samples=args.samples, seed=args.seed)
        expect = real.spec["regulation"].get("expect")
        res["expected_class"] = expect
        res["matches_expectation"] = expect is None or res["class"] == expect or (
            isinstance(expect, list) and res["class"] in expect)
        if not res["matches_expectation"]:
            bad.append(f.name)
        rows.append({"spec": str(f), "report": _report(res, real), "spec_regulation": real.spec["regulation"]})
        print(("   " if res["matches_expectation"] else "!! ") + summary(res), flush=True)
    table = _survey_tables(rows)
    write_report(Path(args.out_dir), "survey", {"command": "regulation survey", "samples": args.samples,
                                               "seed": args.seed, "not_as_expected": bad,
                                               "results": [r["report"] for r in rows]}, "\n".join(table) + "\n")
    print(f"{len(rows)} specifications, {len(bad)} not as expected")
    return 1 if bad else 0


def cmd_check(args) -> int:
    from .certificate import check
    from .spec import load
    report = json.loads(Path(args.report).read_text(encoding="utf-8"))
    spec_path = args.spec or report["card"].get("spec")
    res = check(report, load(spec_path), seed=args.seed)
    for c in res["checks"]:
        print(("passed  " if c["passed"] else "FAILED  ") + c["check"])
    print("certificate " + ("passed" if res["passed"] else "failed"))
    return 0 if res["passed"] else 1


def add_parser(sub) -> None:
    reg = sub.add_parser("regulation", help="The return of an output to its set point after a step of an input: "
                                            "integrators, static gains and step responses.")
    rsub = reg.add_subparsers(dest="regulation_command", required=True)
    p = rsub.add_parser("card", help="Classify one specification and write its report and certificate.")
    p.add_argument("spec")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--samples", type=int, default=32, help="parameter points for the robustness of the class")
    p.add_argument("--seed", type=int, default=0)
    p.set_defaults(func=cmd_card)
    p = rsub.add_parser("survey", help="Classify every specification in the given files or folders and compare with "
                                       "the class each one expects.")
    p.add_argument("paths", nargs="+")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--samples", type=int, default=32)
    p.add_argument("--seed", type=int, default=0)
    p.set_defaults(func=cmd_survey)
    p = rsub.add_parser("check", help="Check the certificate of a report against its specification, without a search.")
    p.add_argument("report")
    p.add_argument("--spec", help="the specification (default: the path recorded in the report)")
    p.add_argument("--seed", type=int, default=1)
    p.set_defaults(func=cmd_check)
