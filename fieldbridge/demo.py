"""``fieldbridge demo``: three calculations from the example specifications, written to one page with their figures.

  1  memory card              the genetic toggle switch: stable states, the write point and its normal form, the
                              writing protocols and the law of loss
  2  phase locking            eight oscillator models from seven fields driven periodically: the derivation in each
                              model as a chain of verified transformations, the step at which a derivation stops, and
                              the invariants of the end points
  3  Bloch rotation           ten quantum carriers (spins, atoms in two wells, exchange chains, Cooper pairs): the
                              rotation derived from each carrier's own Hamiltonian, or the term that obstructs it
  4  adding a material        the commands for a new material specification (CONTRIBUTING.md)

Everything is calculated on this machine from the example specifications. The page (index.html, self-contained,
with the figures embedded) and the report of every step are written to --out-dir; demo.md repeats the page in
Markdown. The calculations need the memory extra: pip install -e '.[memory]'.
"""
from __future__ import annotations

import html
import json
import time
from pathlib import Path
from typing import Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
BOUNDARY = ("Every statement on this page is calculated from a supplied model. It does not establish that the model "
            "describes a physical system, and it does not establish novelty; both need comparison with experiment "
            "and with the literature.")


def _material(out: Path, seed: int) -> Dict:
    from .memory import discovery, predict, spec
    from .memory.visual import card_figure
    real = spec.load(EXAMPLES / "memory" / "toggle.json")
    pred = predict.predict(real, np.random.default_rng(seed))
    card = discovery.evaluate(real, np.random.default_rng(seed + 1), quick=True)
    card["id"] = "toggle"
    comparison = predict.compare(pred, card)
    step = out / "1_material"
    step.mkdir(parents=True, exist_ok=True)
    (step / "card.json").write_text(json.dumps({"structure": pred, "card": card, "comparison": comparison},
                                               indent=1, default=float) + "\n", encoding="utf-8")
    card_figure(card, step / "card.png")
    event = next((e for e in card["construct"]["events"] if e["source"] == "control"), {})
    return {"figure": step / "card.png", "verdict": card["verdict"], "states": card["states"]["count"],
            "write_point": event.get("value"), "write_kind": str(event.get("kind", "")).split(":")[0],
            "threshold": card.get("threshold"), "agreement": comparison.get("all_consistent"),
            "prediction": pred["predictions"].get("write", {}).get("prediction", "")}


def _mechanism(out: Path, seed: int, law: bool) -> Dict:
    from .memory import codiscovery, spec
    from .memory.visual import codiscovery_figure
    reals = [spec.load(p) for p in sorted((EXAMPLES / "memory" / "oscillators").glob("*.json"))]
    report = codiscovery.codiscover(reals, np.random.default_rng(seed), target="phase-locking", check_law=law)
    step = out / "2_mechanism"
    step.mkdir(parents=True, exist_ok=True)
    (step / "codiscover.json").write_text(json.dumps(report, indent=1, default=float) + "\n", encoding="utf-8")
    codiscovery_figure(report, step / "codiscover.png")
    rows = [{"name": r.get("name"), "field": r.get("field"), "status": r.get("status"), "word": r.get("word"),
             "obstruction": r.get("obstruction_short") or r.get("obstruction")} for r in report["rows"]]
    return {"figure": step / "codiscover.png", "summary": report["summary"], "rows": rows}


def _quantum(out: Path) -> Dict:
    from .quantum import language as ql
    reals = [ql.load(p) for p in sorted((EXAMPLES / "quantum").glob("*.json"))]
    report = ql.codiscover(reals)
    step = out / "3_quantum"
    step.mkdir(parents=True, exist_ok=True)
    (step / "codiscover.md").write_text(ql.markdown(report), encoding="utf-8")
    figure = None
    try:
        from .quantum.figure import codiscovery_figure
        codiscovery_figure(report, step / "codiscover.png")
        figure = step / "codiscover.png"
    except ImportError:  # matplotlib is missing
        pass
    return {"figure": figure, "summary": report["summary"]}


def _img(path) -> str:
    from .memory.visual import _img64
    return _img64(path) if path and Path(path).exists() else ""


def _page(results: Dict, out: Path) -> None:
    from .memory.visual import CSS, FONTS
    esc = html.escape
    m, c, q = results["material"], results["mechanism"], results["quantum"]
    s = c["summary"]
    reached = [r for r in c["rows"] if str(r["status"]).startswith("reached")]
    stopped = [r for r in c["rows"] if not str(r["status"]).startswith("reached")]
    qs = q["summary"]
    parts = ["<main>",
             "<header style='display:grid;gap:12px'><h1>FieldBridge demonstration</h1>"
             "<p class='lead'>Three calculations from the example specifications of the repository: the memory card of "
             "the genetic toggle switch, phase locking derived in eight oscillator models from seven fields, and the "
             "Bloch rotation on ten quantum carriers. The last section lists the commands for adding a material."
             "</p></header>",
             "<section class='card'><div class='card-head'><h2>1. Memory card of the genetic toggle switch</h2>"
             "<span class='src'>Gardner, Cantor and Collins (2000)</span></div>"
             "<p>Two genes repress each other. Before any simulation, the exchange of the two genes predicts that if "
             "the symmetric state loses stability along the mode that the exchange reverses, it does so at a "
             "pitchfork, where a weak bias decides which gene is written. The calculation finds "
             f"<b>{m['states']}</b> stable states and a <b>{esc(m['write_kind'])}</b> at promoter strength "
             f"<b>{m['write_point']:.2f}</b>; a uniform field removes the stored state at <b>{m['threshold']:.3f}</b>. "
             f"Every structural prediction agrees with the calculation: <b>{'yes' if m['agreement'] else 'no'}</b>.</p>"
             f"<div class='plate'><img alt='Memory card of the genetic toggle switch' src='{_img(m['figure'])}'></div>"
             "<p class='lead'>(a) flow and stable states, (b) states against the promoter strength, with the write "
             "point, (c) writing protocols and retention without a field, (d) the drift along the unstable direction "
             "at the write point, a pitchfork. Module 1 of the tutorial explains each panel.</p></section>",
             "<section class='card'><div class='card-head'><h2>2. Phase locking in eight oscillator models</h2>"
             "</div><p>Phase locking is derived in oscillators from electronics, chemistry, neuroscience, ecology, "
             "chronobiology, superconductivity and computing hardware. Each derivation is a chain of verified "
             "transformations (a stable cycle, its phase response, the averaged drive) that ends on the Adler "
             "equation d&psi;/dt = &Delta;&omega; &minus; K sin &phi;. "
             f"<b>{len(reached)}</b> of <b>{len(c['rows'])}</b> models reach it"
             + (f"; in the others the derivation stops at a stated step ({esc('; '.join(sorted({str(r['obstruction'])[:60] for r in stopped})))})"
                if stopped else "")
             + ". The ratio at which an oscillator locks is set by a symmetry of the model that the drive respects.</p>"
             f"<div class='plate'><img alt='Phase locking derived in oscillators from different fields' "
             f"src='{_img(c['figure'])}'></div>"
             "<p class='lead'>Module 9 of the tutorial explains the derivations, the obstructions and the invariants."
             "</p></section>",
             "<section class='card'><div class='card-head'><h2>3. The Bloch rotation on ten quantum carriers</h2></div>"
             "<p>The Bloch rotation of a spin is detached from its carrier: the Hamiltonian and the observable generate "
             "the Lie algebra su(2), so the signal follows the Rabi law on every carrier where that closure holds. "
             f"The derivation reaches it on <b>{qs.get('reached')}</b> of <b>{qs.get('reached', 0) + qs.get('obstructed', 0)}</b> "
             "carriers, from their own Hamiltonians; where it stops, the term that enlarges the algebra is named.</p>"
             + (f"<div class='plate'><img alt='The Bloch rotation on different carriers' src='{_img(q['figure'])}'></div>"
                if q["figure"] else "")
             + "<p class='lead'>Chapter 24 of the tutorial shows the language on these spins.</p></section>",
             "<section class='card'><div class='card-head'><h2>4. Adding a material specification</h2></div>"
             "<p>A material is one JSON file: its variables, equations, parameters, the control an experiment varies, "
             "the closure, the observable, the assumptions and the source of the equations.</p>"
             "<pre class='mono' style='white-space:pre-wrap;background:var(--panel);border:1px solid var(--rule);"
             "border-radius:8px;padding:12px'>python3 -B -m fieldbridge memory new my_material --carrier orthant\n"
             "python3 -B -m fieldbridge memory check examples/memory/my_material.json\n"
             "python3 -B -m fieldbridge memory card examples/memory/my_material.json --out-dir build/my_material</pre>"
             "<p>CONTRIBUTING.md describes each field and how to open a pull request; docs/materials.md lists the "
             "materials already described.</p></section>",
             f"<p class='lead'>{esc(BOUNDARY)}</p>", "</main>"]
    head = f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'><title>FieldBridge demonstration</title>{FONTS}<style>{CSS}</style></head><body>"
    (out / "index.html").write_text(head + "".join(parts) + "</body></html>", encoding="utf-8")
    md = ["# FieldBridge demonstration", "", "`index.html` contains the same results with the figures embedded.", "",
          "## 1. Memory card of the genetic toggle switch", "",
          f"Genetic toggle switch: {m['states']} stable states; {m['write_kind']} at promoter strength "
          f"{m['write_point']:.2f}; a uniform field removes the stored state at {m['threshold']:.3f}. Structural "
          f"predictions agree with the calculation: {'yes' if m['agreement'] else 'no'}.", "",
          "![card](1_material/card.png)", "", "## 2. Phase locking in eight oscillator models", "",
          f"Phase locking reached in {len(reached)} of {len(c['rows'])} oscillator models.", "",
          "| Model | Field | Derivation | Where it stops |", "| --- | --- | --- | --- |"]
    md += [f"| {r['name']} | {r['field']} | {r['word'] or ''} | {r['obstruction'] or ''} |" for r in c["rows"]]
    md += ["", "![phase locking](2_mechanism/codiscover.png)", "", "## 3. The Bloch rotation on ten quantum carriers", "",
           f"The Bloch rotation is reached on {qs.get('reached')} carriers and obstructed on {qs.get('obstructed')}.", ""]
    if q["figure"]:
        md += ["![Bloch rotation](3_quantum/codiscover.png)", ""]
    md += ["## 4. Adding a material specification", "", "CONTRIBUTING.md describes the specification fields and the "
           "checks.", "", BOUNDARY, ""]
    (out / "demo.md").write_text("\n".join(md), encoding="utf-8")


def run(out_dir, seed: int = 20260923, law: bool = False, log=print) -> Dict:
    if not (EXAMPLES / "memory" / "toggle.json").exists():
        raise SystemExit("The demonstration reads the example specifications of a clone of the repository; install "
                         "it from the clone with: pip install -e '.[memory]'")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    results, t0 = {}, time.time()
    log("1/3 memory card of the genetic toggle switch ...")
    results["material"] = _material(out, seed)
    log(f"    {time.time() - t0:.0f} s. 2/3 phase locking in oscillators from different fields ...")
    results["mechanism"] = _mechanism(out, seed, law)
    log(f"    {time.time() - t0:.0f} s. 3/3 the Bloch rotation on quantum carriers ...")
    results["quantum"] = _quantum(out)
    _page(results, out)
    log(f"    {time.time() - t0:.0f} s. Open {out / 'index.html'}")
    return results


def add_parser(sub) -> None:
    p = sub.add_parser("demo", help="Three calculations from the example specifications, written to one page with "
                                    "their figures.")
    p.add_argument("--out-dir", default="build/demo")
    p.add_argument("--seed", type=int, default=20260923)
    p.add_argument("--law", action="store_true",
                   help="Also measure the locking law of every oscillator (slower).")
    p.set_defaults(func=lambda args: (run(args.out_dir, seed=args.seed, law=args.law), 0)[1])
