"""Collect the hardware survey (survey.py) into examples/hardware/REPORT.md: the coverage of the census, the three
derivations of every model, the laws at their obstructions and the single changes that reach a write.

    python3 -B examples/hardware/report.py [--results examples/hardware/results]
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPECTED = {"symmetric-write": 1.3313353638003895, "threshold-write": 1.0187929, "phase-locking": 1.0}
NAMES = {"symmetric-write": "symmetric write", "threshold-write": "threshold write", "phase-locking": "phase locking"}


def derivation_cell(t: str, d: dict) -> str:
    status = str(d.get("status") or "")
    lc = d.get("law_constant") or {}
    if status.startswith("reached"):
        c, e = lc.get("constant"), lc.get("stderr")
        if c is None:
            return "reached, no law"
        known = abs(c - EXPECTED[t]) <= max(3 * (e or 0.0), 0.03 * EXPECTED[t])
        return f"{'known law' if known else 'law differs'}, {c:.3f}" + (f" ± {e:.3f}" if e else "")
    return str(d.get("obstruction_short") or d.get("obstruction") or "")[:70]


def law_rows(comp: dict):
    """One row per law derived at an obstruction: completion, second parameter, certification, response and offset
    (law and full dynamics) or, for several critical directions, the number of states and the largest deviation."""
    for e in comp.get("write_points", []):
        laws = [("", e["law"])] if "law" in e else [(u["second"], u.get("law", {"status": "no cusp found"}))
                                                     for u in e.get("unfoldings", [])]
        if not laws:
            yield {"completion": e["completion"], "second": "", "status": "no second parameter"}
        for second, L in laws:
            row = {"completion": e["completion"], "second": second, "status": str(L.get("status", ""))[:80],
                   "certified": L.get("certified")}
            od = L.get("one_direction") or {}
            for src in ("law", "certificate"):
                xs = od.get(src, [])
                mid = next((x for x in xs if x["bias_floors"] == 1.0), xs[0] if xs else None)
                if mid:
                    row[src] = (round(mid["response_constant"], 2), round(mid["offset"], 2))
            if (L.get("critical_dimension") or 1) > 1 and L.get("rows"):
                row["states"] = len(L.get("stable_states", []))
                row["deviation"] = round(max(r.get("deviation", 0) for r in L["rows"]), 3)
            yield row


def law_text(r: dict) -> str:
    head = r["completion"] + (f" ({r['second']})" if r.get("second") else "")
    if not r["status"].startswith("derived"):
        return f"{head}: {r['status']}"
    cert = "certified" if r.get("certified") else "not certified"
    if "states" in r:
        return f"{head}: {r['states']} states, largest deviation {r['deviation']}, {cert}"
    if "law" in r:
        return (f"{head}: response {r['law'][0]}, offset {r['law'][1]}; full dynamics {r['certificate'][0]}, "
                f"{r['certificate'][1]}; {cert}" if "certificate" in r else f"{head}: response {r['law'][0]}, "
                f"offset {r['law'][1]}; {cert}")
    return f"{head}: {cert}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=HERE / "results")
    args = ap.parse_args()
    census = json.loads((HERE / "census.json").read_text(encoding="utf-8"))
    res = {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in sorted(args.results.glob("*.json"))
           if not p.name.endswith(".partial.json")}
    cov = Counter(e["classical_write"] for e in census["census"])
    lines = ["# Mechanisms of quantum-computing hardware", "",
             "Written by `python3 -B examples/hardware/report.py` from the results of `survey.py`.", "",
             f"The census lists {len(census['census'])} hardware elements described in reviews of superconducting, "
             f"semiconductor, trapped-ion, neutral-atom, photonic and optomechanical platforms: {cov['yes']} have a "
             f"classical write, {cov['partly']} have one in part and {cov['no']} have none (gates, codes and linear "
             f"amplifiers store or process information without a write). {len(res)} models represent the elements "
             "with a write; the status of every element is in `census.json`.", "",
             "## Derivations and laws at the obstructions", "",
             "| model | symmetric write | threshold write | phase locking | laws at the obstructions | a single change "
             "that reaches a write |", "| --- | --- | --- | --- | --- | --- |"]
    tally = Counter()
    offsets, planes = [], []
    for m, r in res.items():
        cells = []
        for t in EXPECTED:
            c = derivation_cell(t, r["derivations"][t])
            cells.append(c)
            tally[f"{NAMES[t]}: " + ("known law" if c.startswith("known law") else "reached" if c.startswith(("reached", "law differs"))
                                     else "stops at an obstruction")] += 1
        rows = list(law_rows(r["completion"]))
        for x in rows:
            tally["law at an obstruction: " + ("certified" if x["status"].startswith("derived") and x.get("certified")
                                              else "not certified" if x["status"].startswith("derived")
                                              else x["status"].split(":")[0][:40])] += 1
            if x["completion"] == "second parameter" and "certificate" in x:
                offsets.append((m, x))
            if "states" in x:
                planes.append((m, x))
        laws = "; ".join(law_text(x) for x in rows) or "no write point along the control"
        steps = [s for s in r["changes"]["singles"] if s["reaches"]]
        ch = "; ".join(f"{s['what']} ({', '.join(NAMES[t] for t in s['reaches'])})" for s in steps[:2]) or ""
        lines.append(f"| {m} | {cells[0]} | {cells[1]} | {cells[2]} | {laws} | {ch} |")
    lines += ["", "## Counts", ""] + [f"- {k}: {v}" for k, v in sorted(tally.items())]
    lines += ["", "## Sweeps through a cusp: response to the bias and intrinsic offset", "",
              "Response constant and offset at a bias of one noise floor, from the law of the completed model and "
              "from the full stochastic dynamics.", "", "| model | second parameter | law | full dynamics |",
              "| --- | --- | --- | --- |"]
    lines += [f"| {m} | {x['second']} | {x['law'][0]}, {x['law'][1]} | {x['certificate'][0]}, {x['certificate'][1]} |"
              for m, x in offsets if "law" in x]
    lines += ["", "## Writes into a critical plane", "", "| model | written states | largest deviation | certified |",
              "| --- | ---: | ---: | --- |"]
    lines += [f"| {m} | {x['states']} | {x['deviation']} | {'yes' if x.get('certified') else 'no'} |" for m, x in planes]
    lines += ["", "Calculated statements about the classical images of the hardware elements, not about the devices; "
              "the quantum dynamics at a few photons can differ (Module 15).", ""]
    (HERE / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines[:12]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
