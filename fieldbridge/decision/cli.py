"""``fieldbridge decision``: predictions, cards and surveys of populations that decide.

    fieldbridge decision predict examples/decision/synchronization/van_der_pol.json
    fieldbridge decision card examples/decision/write/honeybees_pais2013.json --simulate --replicas 2000
    fieldbridge decision survey examples/decision --out-dir out/decision
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

from ..core.provenance import provenance, write_report

BOUNDARY = ("Calculated statements about a supplied model of a population. Where the source has no noise, a ramp or a "
            "bias, these are assumptions stated in the specification, and the probabilities and times depend on them. "
            "The stochastic checks run the mean-field equations with a noise that falls as 1/N (or the units "
            "themselves); they are not the individual-based models of the sources. The statements do not establish "
            "that the model describes a physical system, and they do not establish novelty; both need comparison "
            "with experiment and with the literature.")
HERE = Path(__file__).resolve().parent
CLASS_TEXT = {"synchronizes": "synchronizes above K_c", "no-onset": "no onset of synchrony",
              "follows-the-bias": "follows the bias (swept write)", "set-by-the-sample": "set by the sample",
              "reflection-seed": "passage from a one-component seed", "rotation-seed": "passage from a two-component seed"}


def _plain(x):
    import numpy as np
    if isinstance(x, dict):
        return {str(k): _plain(v) for k, v in x.items() if not str(k).startswith("_") and not callable(v)}
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        return x.item()
    if isinstance(x, np.ndarray):
        return _plain(x.tolist())
    return x


def _f(v, fmt=".4g"):
    return "—" if v is None else format(v, fmt)


def markdown(res: Dict) -> str:
    pr = res["prediction"]
    red = pr["reduction"]
    lines = [f"# Decision card: {res['name']}", "", f"Field: {res['field'] or 'none'}. Source: {res['source']}", "",
             f"Kind: {res['kind']}; target: {res['target']}.", "", f"## Class: {CLASS_TEXT[res['class']]}", ""]
    if res["target"] == "synchronization":
        lines += [f"- Unit: frequency {red['omega']:.5g}, Floquet rate kappa = {red['kappa']:.4g}; coupling function "
                  f"a1 = {red['a1']:+.4g}, b1 = {red['b1']:+.4g}.",
                  f"- Frequency spread from the parameter spread {red['param_spread']:.4g} (d omega/dp = "
                  f"{red['domega_dp']:.4g}).",
                  f"- K_c = {_f(red['K_c'])}, Omega = {_f(red['Omega_c'])}.", ""]
        if red["rows"]:
            lines += ["| K/K_c | K | mu | K/kappa |", "| --- | --- | --- | --- |"]
            lines += [f"| {r['factor']:g} | {r['K']:.4g} | {r['mu']:.4g} | {r['K_over_kappa']:.3g} |" for r in red["rows"]]
            lines.append("")
    else:
        keys = ("c_star", "a", "b", "h_s")
        lines += ["- Reduction at the threshold: " + ", ".join(f"{k} = {red[k]:.5g}" for k in keys if k in red) + ".", ""]
    if res["target"] == "collective-write":
        rows = pr["rows"]
        if "s_frozen" in rows[0]:
            lines += ["| N | r | z | h | P | sigma_thermal | s_frozen |", "| --- | --- | --- | --- | --- | --- | --- |"]
            lines += [f"| {r['N']:g} | {r['r']:.4g} | {r['z']:g} | {r['h']:.3e} | {r['P']:.4f} | "
                      f"{r['sigma_thermal']:.3e} | {r['s_frozen']:.3e} |" for r in rows]
        else:
            lines += ["| N | r | z | h | P (passage) | P (closed form) | Λ | window |",
                      "| --- | --- | --- | --- | --- | --- | --- | --- |"]
            lines += [f"| {r['N']:g} | {r['r']:.4g} | {r['z']:g} | {r['h']:.3e} | {r['P']:.4f} | {r['P_linear']:.4f} | "
                      f"{r['Lambda']:.3g} | {r['window']:.3g} |" for r in rows]
        lines.append("")
    if res["target"] == "seeded-passage":
        lines += [f"- Leading eigenspace after the step: d = {pr['d']}; rate at the end {pr['rate_end']:.4g}.",
                  f"- Window law: Lambda(tau_90) - Lambda(tau_10) = {pr['window_law']:.4f}"
                  + (f" (constant rate: {pr['window_time_at_end']:.4g} in time)." if pr['window_time_at_end'] else "."),
                  ""]
    if "simulation" in res:
        sim = res["simulation"]
        lines += ["## Stochastic check", ""]
        if res["target"] == "synchronization":
            lines.append(f"{sim['units']} full units at K = {sim['factor']:g} K_c: growth rate of synchrony "
                         f"{sim['mu_measured']:.4g} against {sim['mu_reduction']:.4g} (ratio {sim['ratio']:.3f}, "
                         f"K/kappa {sim['K_over_kappa']:.2g}; the excess is of order K/kappa).")
        elif res["target"] == "collective-write":
            lines += ["| N | r | z | measured | law | within 4 stderr + 0.02 |", "| --- | --- | --- | --- | --- | --- |"]
            lines += [f"| {r['N']:g} | {r['r']:.4g} | {r['z']:g} | {r['measured']:.4f} ± {r['stderr']:.4f} | "
                      f"{r['P']:.4f} | {'yes' if r['agrees'] else 'no'} |" for r in sim["rows"]]
        else:
            lines += ["| N | window in Lambda | law | ratio | Lambda(tau_50) | within gate |",
                      "| --- | --- | --- | --- | --- | --- |"]
            for r in sim["rows"]:
                if "error" in r:
                    lines.append(f"| {r['N']:g} | {r['error']} | | | | |")
                else:
                    lines.append(f"| {r['N']:g} | {r['Lambda_window']:.4f} ± {r['Lambda_window_se']:.4f} | "
                                 f"{r['law']:.4f} | {r['ratio']:.4f} | {r['Lambda_median']:.3f} | "
                                 f"{'yes' if r['agrees'] else 'no'} |")
            if "median_slope" in sim:
                m = sim["median_slope"]
                lines += ["", f"Lambda(tau_50) per ln N: {m['slope']:.4f} ± {m['slope_se']:.4f} (law 1/2)."]
        lines.append("")
    lines += ["## Boundary", "", BOUNDARY, ""]
    return "\n".join(lines)


def _load(path):
    from .spec import load
    return load(path)


def cmd_predict(args) -> int:
    from .predict import predict
    pop = _load(args.spec)
    pr = predict(pop)
    red = pr["reduction"]
    if pop.target == "synchronization":
        extra = f"K_c = {_f(red['K_c'])}, a1 = {red['a1']:+.3g}, b1 = {red['b1']:+.3g}"
    elif pop.target == "collective-write":
        extra = f"c* = {red['c_star']:.5g}, first condition P = {pr['rows'][0]['P']:.3f}"
    else:
        extra = f"d = {pr['d']}, rate at the end {pr['rate_end']:.4g}, window law {pr['window_law']:.4f}"
    print(f"{pop.name}: {CLASS_TEXT[pr['class']]} ({extra})")
    return 0


def cmd_card(args) -> int:
    from .card import card
    pop = _load(args.spec)
    res = card(pop, simulate=args.simulate, replicas=args.replicas, seed=args.seed)
    report = _plain(dict(res, provenance=provenance(HERE, BOUNDARY, pop.spec)))
    name = Path(args.spec).stem
    write_report(Path(args.out_dir), name, report, markdown(res))
    print(f"{pop.name}: {CLASS_TEXT[res['class']]}; report in {Path(args.out_dir) / (name + '.md')}")
    return 0


def survey_rows(folder: Path) -> List[Dict]:
    from .predict import predict
    rows = []
    for path in sorted(Path(folder).rglob("*.json")):
        pop = _load(path)
        pr = predict(pop)
        rows.append({"path": str(path), "name": pop.name, "field": pop.field, "target": pop.target,
                     "class": pr["class"], "expect": pop.expect})
    return rows


def cmd_survey(args) -> int:
    rows = survey_rows(Path(args.folder))
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "survey.json").write_text(json.dumps(_plain({"rows": rows}), indent=2) + "\n", encoding="utf-8")
    md = ["| Population | Field | Target | Class |", "| --- | --- | --- | --- |"]
    md += [f"| {r['name']} | {r['field'] or '—'} | {r['target']} | {CLASS_TEXT[r['class']]} |" for r in rows]
    text = "\n".join(md) + "\n"
    (out / "survey.md").write_text(text, encoding="utf-8")
    print(text)
    bad = [r["path"] for r in rows if r["expect"] is not None and r["expect"] != r["class"]]
    if bad:
        print("Classes differ from the expected ones: " + ", ".join(bad))
        return 1
    return 0


def add_parser(sub) -> None:
    p = sub.add_parser("decision", help="Collective decisions: synchronization, the swept write, passage from a seed")
    ds = p.add_subparsers(dest="decision_command", required=True)
    pp = ds.add_parser("predict", help="Reduction, law and class from the equations")
    pp.add_argument("spec")
    pp.set_defaults(func=cmd_predict)
    pc = ds.add_parser("card", help="The decision card, with a stochastic check on request")
    pc.add_argument("spec")
    pc.add_argument("--simulate", action="store_true", help="Run the stochastic check of the target")
    pc.add_argument("--replicas", type=int, default=2000)
    pc.add_argument("--seed", type=int, default=1)
    pc.add_argument("--out-dir", default="build/decision")
    pc.set_defaults(func=cmd_card)
    ps = ds.add_parser("survey", help="Predictions for every specification in a folder")
    ps.add_argument("folder")
    ps.add_argument("--out-dir", default="build/decision/survey")
    ps.set_defaults(func=cmd_survey)
