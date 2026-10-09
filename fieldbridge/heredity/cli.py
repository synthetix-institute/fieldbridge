"""``fieldbridge heredity``: predictions, cards and surveys of bodies that grow and divide.

    fieldbridge heredity predict examples/heredity/turing_painter1999.json
    fieldbridge heredity card examples/heredity/filament_baczynski2007.json --simulate --lineages 1000 --out-dir out/fil
    fieldbridge heredity survey examples/heredity --out-dir out/heredity
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

from ..core.provenance import provenance, write_report

BOUNDARY = ("Calculated statements about a supplied model of a body that grows and divides. Growth, division and, "
            "where the source has none, the noise are assumptions stated in the specification; the inherited fraction "
            "depends on them. The statements do not establish that the model describes a physical system, and they do "
            "not establish novelty; both need comparison with experiment and with the literature.")
HERE = Path(__file__).resolve().parent
CLASS_TEXT = {"inherited-through-threshold": "inherited through a threshold",
              "kept-above-threshold": "kept above the threshold",
              "lost-in-the-dip": "lost in the dip",
              "threshold-moved": "threshold moved by division"}


def _plain(x):
    import numpy as np
    if isinstance(x, dict):
        return {str(k): _plain(v) for k, v in x.items() if not str(k).startswith("_")}
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        return x.item()
    if isinstance(x, np.ndarray):
        return _plain(x.tolist())
    return x


def markdown(res: Dict) -> str:
    pr = res["prediction"]
    g = res["growth"]
    lines = [f"# Heredity card: {res['name']}", "", f"Field: {res['field'] or 'none (written without a field)'}. "
             f"Source: {res['source']}", "",
             f"Body: {res['body']}; size `{res['size']}`, {g['law']} growth at the rate {g['rate']:g}; division at "
             f"{pr['L_div']:.4g}; noise {res['noise']}.", "",
             f"## Class: {CLASS_TEXT[res['class']]}", "",
             f"- Threshold of the symmetric state: L_c = {pr['L_c']:.5g}; reduction a = {pr['a']:.4g}, b = {pr['b']:+.4g} "
             f"({'supercritical' if pr['supercritical'] else 'subcritical: the law does not apply'}).",
             f"- Dip: {'the daughter is born at ' + format(pr['L_div'] / 2, '.4g') + ' < L_c' if pr['dip'] else 'none, the daughter is born above L_c'}.",
             (f"- Linear gain over a generation ln G = {pr['lnG']:.3g} "
              f"({'the order survives a generation without noise' if pr['lnG'] > 0 else 'the order is lost without noise'})."
              if pr["lnG"] is not None else "- Linear gain over a generation: not needed without a dip."),
             f"- Order at successive divisions of a lineage without noise: "
             f"{', '.join(f'{v:.3g}' for v in pr['order_at_divisions'])}; displacement of the daughter from the "
             f"symmetric state outside the order: {pr['displacement']:.2g}.",
             f"- Lambda = a ramp / (b D_s) = {pr['Lambda']:.3g} (the law needs Lambda >> 1).", ""]
    if "law" in pr:
        law = pr["law"]
        lines += ["## Law for a daughter", "",
                  f"Order at birth {law['order_at_birth']:.4g}, at the crossing {law['order_at_crossing']:.4g} (from the "
                  f"body's own equations); noise amplitude sigma_c = {law['sigma_linear']:.4g} (linear passage), "
                  f"{law['sigma_history']:.4g} (eigenvalue history).", "",
                  f"P = Phi(phi_c / sigma_c) = {law['P_body']:.4f} (history {law['P_history']:.4f}; normal form with a "
                  f"linear ramp {law['P_normal_form']:.4f}).", ""]
    if "lineages" in res:
        r = res["lineages"]
        lines += ["## Lineages of the full body", "",
                  f"{r['lineages']} lineages, {r['transfers']} divisions counted: measured {r['measured']:.4f} "
                  f"± {r['stderr']:.4f}; law (body) {r['law_body']:.4f}, history {r['law_history']:.4f}, normal form "
                  f"{r['law_normal_form']:.4f}; within |P - P_law| < 4 stderr + 0.02: "
                  f"{'yes' if r['agrees_body'] else 'no'} (body), {'yes' if r['agrees_history'] else 'no'} (history).", ""]
    lines += ["## Boundary", "", BOUNDARY, ""]
    return "\n".join(lines)


def _load(path):
    from .spec import load
    return load(path)


def cmd_predict(args) -> int:
    from .predict import predict, reduction
    lin = _load(args.spec)
    pr = predict(lin, reduction(lin))
    print(f"{lin.name}: {CLASS_TEXT[pr['class']]} (L_c = {pr['L_c']:.5g}, b = {pr['b']:+.3g}, "
          f"L_div = {pr['L_div']:.4g}, ln G = {pr['lnG'] if pr['lnG'] is None else format(pr['lnG'], '.3g')}, "
          f"Lambda = {pr['Lambda']:.3g}" + (f", P = {pr['law']['P_body']:.3f}" if "law" in pr else "") + ")")
    return 0


def cmd_card(args) -> int:
    from .card import card
    lin = _load(args.spec)
    res = card(lin, simulate=args.simulate, lineages=args.lineages, generations=args.generations, seed=args.seed,
               noise=args.noise_scale)
    report = _plain(dict(res, provenance=provenance(HERE, BOUNDARY, lin.spec)))
    name = Path(args.spec).stem
    write_report(Path(args.out_dir), name, report, markdown(res))
    print(f"{lin.name}: {CLASS_TEXT[res['class']]}; report in {Path(args.out_dir) / (name + '.md')}")
    return 0


def survey_rows(folder: Path) -> List[Dict]:
    from .predict import predict, reduction
    rows = []
    for path in sorted(Path(folder).rglob("*.json")):
        lin = _load(path)
        pr = predict(lin, reduction(lin))
        rows.append({"path": str(path), "name": lin.name, "field": lin.field, "class": pr["class"],
                     "expect": lin.expect, "L_c": pr["L_c"], "b": pr["b"], "L_div": pr["L_div"], "lnG": pr["lnG"],
                     "Lambda": pr["Lambda"], "P_body": pr.get("law", {}).get("P_body"),
                     "P_normal_form": pr.get("law", {}).get("P_normal_form")})
    return rows


def survey_table(rows: List[Dict]) -> str:
    md = ["| Body | Field | Class | L_c | L_div / L_c | ln G | Λ | P (body) | P (normal form) |",
          "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    f = lambda v, fmt: "—" if v is None else format(v, fmt)   # noqa: E731
    for r in rows:
        md.append(f"| {r['name']} | {r['field'] or '—'} | {CLASS_TEXT[r['class']]} | {r['L_c']:.4g} | "
                  f"{r['L_div'] / r['L_c']:.3g} | {f(r['lnG'], '.3g')} | {f(r['Lambda'], '.3g')} | "
                  f"{f(r['P_body'], '.3f')} | {f(r['P_normal_form'], '.3f')} |")
    return "\n".join(md) + "\n"


def cmd_survey(args) -> int:
    rows = survey_rows(Path(args.folder))
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "survey.json").write_text(json.dumps(_plain({"rows": rows}), indent=2) + "\n", encoding="utf-8")
    text = survey_table(rows)
    (out / "survey.md").write_text(text, encoding="utf-8")
    print(text)
    bad = [r["path"] for r in rows if r["expect"] is not None and r["expect"] != r["class"]]
    if bad:
        print("Classes differ from the expected ones: " + ", ".join(bad))
        return 1
    return 0


def add_parser(sub) -> None:
    p = sub.add_parser("heredity", help="Inheritance through growth and division: threshold, law and lineages")
    hs = p.add_subparsers(dest="heredity_command", required=True)
    pp = hs.add_parser("predict", help="Threshold, conditions and class from the equations")
    pp.add_argument("spec")
    pp.set_defaults(func=cmd_predict)
    pc = hs.add_parser("card", help="The heredity card, with lineages of the full body on request")
    pc.add_argument("spec")
    pc.add_argument("--simulate", action="store_true", help="Run lineages of the full body with noise")
    pc.add_argument("--lineages", type=int, default=1000)
    pc.add_argument("--generations", type=int, default=3)
    pc.add_argument("--seed", type=int, default=1)
    pc.add_argument("--noise-scale", type=float, default=1.0, help="Multiply the noise of the specification")
    pc.add_argument("--out-dir", default="build/heredity")
    pc.set_defaults(func=cmd_card)
    ps = hs.add_parser("survey", help="Predictions for every specification in a folder")
    ps.add_argument("folder")
    ps.add_argument("--out-dir", default="build/heredity/survey")
    ps.set_defaults(func=cmd_survey)
