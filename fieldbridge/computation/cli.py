"""``fieldbridge computation``: predictions, cards and surveys of driven bodies.

    fieldbridge computation predict examples/computation/chemotaxis_tu2008.json
    fieldbridge computation card examples/computation/spin_torque_furuta2018.json --noise 1e-3 --out-dir out/sto
    fieldbridge computation card examples/computation/hodgkin_huxley.json --simulate --amplitudes 1,4 --out-dir out/hh
    fieldbridge computation survey examples/computation --out-dir out/computation
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

from ..core.provenance import provenance, write_report

BOUNDARY = ("Calculated statements about a supplied model driven by independent inputs. Capacities depend on the input "
            "law, its amplitude, the hold and the measurement noise, which the report states. They do not establish "
            "that the model describes a physical system, and they do not establish novelty; both need comparison "
            "with experiment and with the literature.")
HERE = Path(__file__).resolve().parent
CLASS_TEXT = {"linear-memory": "linear memory only", "odd-capacity": "odd degrees only",
              "nonlinear-capacity": "nonlinear capacity", "integrating": "no fading memory"}


def _plain(x):
    import numpy as np
    if isinstance(x, dict):
        return {str(k): _plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        return x.item()
    if isinstance(x, np.ndarray):
        return _plain(x.tolist())
    return x


def _degrees(d: Dict) -> str:
    return ", ".join(f"degree {k}: {v:.4g}" for k, v in sorted(d.items(), key=lambda kv: int(kv[0])))


def markdown(res: Dict) -> str:
    pr = res["prediction"]
    lines = [f"# Computation card: {res['name']}", "", f"Field: {res['field']}. Source: {res['source']}", "",
             f"Input `{res['input']}`, {res['law']} with amplitude {res['amplitude']:g} around {res['offset']:g}, held "
             f"for {res['hold']:g}; measured: {', '.join(res['observables'])}"
             + (f", {res['virtual_nodes']} times per interval" if res["virtual_nodes"] > 1 else "") + ".", "",
             f"## Class: {CLASS_TEXT[res['class']]}", "",
             f"- Memory: {pr['class']}; slowest rate {pr['slowest_rate']:.4g} per unit time, over the "
             f"{pr['modes_reached_and_seen']} of {pr['modes']} modes that the input reaches and the observables see.",
             f"- Rank of the linear response n_lin = {pr['n_lin']} of {pr['signals']} measured signals (relative "
             f"singular values {', '.join(f'{v:.2g}' for v in pr['singular_values'][:8])}).",
             f"- Degree-1 profile at small amplitude, delays 0-9: "
             f"{', '.join(f'{v:.3f}' for v in pr['profile_degree_1'][:10])}.",
             f"- Odd about the steady state: {'yes (even degrees vanish)' if pr['odd'] else 'no'}; linear: "
             f"{'yes (degree 1 only)' if pr['linear'] else 'no'}.", ""]
    for run in res["runs"]:
        lines.append(f"## Amplitude {run['amplitude']:g}" + (f", measurement noise {res['noise']:g}" if res["noise"] else ""))
        lines.append("")
        if "exact" in run:
            lines.append(f"- Exact (grid {run['exact']['grid']}): {_degrees(run['exact']['by_degree'])}; total "
                         f"{run['exact']['total']:.4g}.")
        if "simulated" in run:
            s = run["simulated"]
            lines.append(f"- Simulated (T = {s['T']}, rank {s['rank']}, threshold {s['threshold']:.2g}): "
                         f"{_degrees(s['by_degree'])}; total {s['total']:.4g}.")
        lines.append("")
    lines += ["## Boundary", "", BOUNDARY, ""]
    return "\n".join(lines)


def _delays(text: str) -> Dict[int, int]:
    return {int(k): int(v) for k, v in (kv.split(":") for kv in text.split(","))}


def cmd_predict(args) -> int:
    from . import spec as cspec
    from .predict import predict
    body = cspec.load(args.spec)
    pr = predict(body)
    print(f"{body.name}: {CLASS_TEXT[pr['structure']]}; memory {pr['class']}, slowest rate {pr['slowest_rate']:.4g}; "
          f"n_lin {pr['n_lin']} of {pr['signals']} signals; modes reached and seen {pr['modes_reached_and_seen']} of "
          f"{pr['modes']}; odd {pr['odd']}; linear {pr['linear']}")
    return 0


def cmd_card(args) -> int:
    from . import spec as cspec
    from .card import card
    body = cspec.load(args.spec)
    amps = [float(a) for a in args.amplitudes.split(",")] if args.amplitudes else None
    res = card(body, delays=_delays(args.delays), amplitudes=amps, simulate_runs=args.simulate, noise=args.noise,
               streams=args.streams, length=args.length, seed=args.seed)
    res = _plain(res)
    print(markdown(res))
    if args.out_dir:
        spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
        res["provenance"] = provenance(HERE, BOUNDARY, spec)
        write_report(Path(args.out_dir), "computation", res, markdown(res))
    return 0


def cmd_survey(args) -> int:
    from . import spec as cspec
    from .card import card
    rows: List[Dict] = []
    paths = sorted(p for p in Path(args.folder).rglob("*.json"))
    for path in paths:
        body = cspec.load(path)
        res = _plain(card(body, delays=_delays(args.delays)))
        pr = res["prediction"]
        ex = res["runs"][0].get("exact")
        rows.append({"path": str(path), "name": res["name"], "field": res["field"], "class": res["class"],
                     "expect": res.get("expect"), "memory": pr["class"], "slowest_rate": pr["slowest_rate"],
                     "modes": f"{pr['modes_reached_and_seen']} of {pr['modes']}", "n_lin": pr["n_lin"],
                     "signals": pr["signals"], "odd": pr["odd"], "linear": pr["linear"],
                     "exact": ex["by_degree"] if ex else None})
        print(f"{path.name}: {CLASS_TEXT[res['class']]}" + (f" (expected {CLASS_TEXT.get(res.get('expect'), res.get('expect'))})"
                                                          if res.get("expect") else ""), flush=True)
    md = survey_table(rows)
    text = md
    print(text)
    mismatched = [r for r in rows if r["expect"] and r["expect"] != r["class"]]
    if args.out_dir:
        report = {"rows": rows, "provenance": provenance(HERE, BOUNDARY)}
        write_report(Path(args.out_dir), "survey", report, "# Computation survey\n\n" + text)
    if mismatched:
        print("classes different from the expected ones: " + ", ".join(r["path"] for r in mismatched))
        return 1
    return 0


def _cap(v: float) -> str:
    return "0" if abs(v) < 1e-9 else f"{v:.3g}"


def survey_table(rows: List[Dict]) -> str:
    """The markdown table of a survey: one row per specification."""
    md = ["| Body | Field | Class | Memory, slowest rate | Modes reached and seen | n_lin / signals | Odd | Exact capacity by degree |",
          "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in rows:
        exact = "—" if not r["exact"] else ", ".join(f"{k}: {_cap(v)}" for k, v in sorted(r["exact"].items(), key=lambda kv: int(kv[0])))
        memory = "fading, " + f"{r['slowest_rate']:.3g}" if r["memory"] == "fading" else "none (integrating)"
        md.append(f"| {r['name']} | {r['field'] or '—'} | {CLASS_TEXT[r['class']]} | {memory} | {r['modes']} | "
                  f"{r['n_lin']} / {r['signals']} | {'yes' if r['odd'] else 'no'} | {exact} |")
    return "\n".join(md) + "\n"


def add_parser(sub) -> None:
    p = sub.add_parser("computation", help="which functions of an input history a driven body represents")
    s = p.add_subparsers(dest="computation_command", required=True)
    a = s.add_parser("predict", help="the predictions from the equations alone (fast)")
    a.add_argument("spec")
    a.set_defaults(func=cmd_predict)
    a = s.add_parser("card", help="predictions, exact capacities (one or two variables) and, on request, simulated ones")
    a.add_argument("spec")
    a.add_argument("--delays", default="1:40,2:8,3:5", help="degree:largest delay, comma separated")
    a.add_argument("--amplitudes", default=None)
    a.add_argument("--noise", type=float, default=0.0, help="relative measurement noise")
    a.add_argument("--simulate", action="store_true", help="also estimate the capacities from a simulation")
    a.add_argument("--streams", type=int, default=16)
    a.add_argument("--length", type=int, default=None)
    a.add_argument("--seed", type=int, default=1)
    a.add_argument("--out-dir", default=None)
    a.set_defaults(func=cmd_card)
    a = s.add_parser("survey", help="predictions and exact capacities of every specification in a folder")
    a.add_argument("folder")
    a.add_argument("--delays", default="1:40,2:8,3:5")
    a.add_argument("--out-dir", default=None)
    a.set_defaults(func=cmd_survey)
