"""Contributing a material: a template specification, the checks it must pass, and the catalog of materials.

  template(name, carrier)   a specification that loads and runs, with placeholders for everything a contributor
                            must state (question, source, field, assumptions, closure, observable)
  check(spec)               required checks (the specification loads; its question, source, field and assumptions
                            are stated rather than left as placeholders) and recommended ones (a citable reference,
                            a control parameter, closure and observable); then the structural predictions and,
                            unless structure_only, a quick memory card compared with them
  catalog(paths)            one row per material: field, question, source and the calculated memory, for
                            docs/materials.md and docs/materials.json
  same_mechanism(...)       the materials of the catalog whose writing mechanism is the one calculated for a new
                            material: what a contribution connects to in other fields
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Dict, Iterable, List, Union

import numpy as np

from . import spec as fspec

CATALOG_JSON = Path(__file__).resolve().parents[2] / "docs" / "materials.json"

PLACEHOLDERS = ("your name", "author, journal", "what the variables are", "for example soft matter",
                "state each modelling choice", "what an experiment measures", "what is specified or discarded",
                "which states does this material", "todo")
REFERENCE = re.compile(r"10\.\d{4,9}/\S+|arxiv:\s*(?:\d{4}\.\d{4,5}|[a-z\-]+/\d{7})|https?://\S+|\b(?:1[89]|20)\d{2}\b",
                       re.IGNORECASE)

_CARRIERS = {
    "euclid": {"carrier": {"kind": "euclid", "name": "what the variables are (an order parameter, a position)",
                           "variables": ["x"], "scale": 1.0},
               "parameters": {"mu": 1.0}, "drift": {"x": "mu*x - x**3"},
               "control": {"name": "mu", "range": [-1.0, 2.0]}},
    "orthant": {"carrier": {"kind": "orthant", "name": "what the variables are (concentrations, populations)",
                            "variables": ["u", "v"], "scale": 4.0},
                "parameters": {"alpha": 4.0, "n": 2.0},
                "drift": {"u": "alpha/(1 + v**n) - u", "v": "alpha/(1 + u**n) - v"},
                "control": {"name": "alpha", "range": [0.5, 6.0]}},
    "torus": {"carrier": {"kind": "torus", "name": "what the variables are (orientations, phases)",
                          "variables": ["theta"], "scale": 1.0, "period": "2*pi"},
              "parameters": {"k": 1.0}, "drift": {"theta": "-k*sin(2*theta)"},
              "control": {"name": "k", "range": [0.1, 2.0]}},
}


def template(name: str, carrier: str = "euclid") -> Dict:
    """A specification that loads and runs. The dynamics is a placeholder of the chosen carrier kind (a pitchfork,
    two mutually repressing genes, an angle with two preferred orientations); every descriptive field is a
    placeholder that check() reports until it is replaced."""
    if carrier not in _CARRIERS:
        raise ValueError(f"carrier must be one of {sorted(_CARRIERS)}")
    body = json.loads(json.dumps(_CARRIERS[carrier]))
    return {
        "schema": fspec.SCHEMA,
        "question": "Which states does this material store, and how is a state written?",
        "provenance": {"source": "Author, Journal volume, page (year), or arXiv:0000.00000",
                       "kind": "contributed example", "contributor": "your name or GitHub handle"},
        "field": "the field of the material, for example soft matter",
        "kind": "equations",
        "name": name,
        **body,
        "material": [],
        "noise": 0.02,
        "dt": 0.01,
        "closure": "what is specified or discarded to obtain closed equations (bath, boundaries, imposed "
                   "conservation)",
        "observable": "what an experiment measures",
        "assumptions": ["state each modelling choice that a reader could dispute"],
    }


def _placeholder(text: str) -> bool:
    low = str(text).lower()
    return any(p in low for p in PLACEHOLDERS)


def _row(name: str, passed: bool, message: str, required: bool) -> Dict:
    return {"check": name, "passed": bool(passed), "required": required, "message": message}


def check(source: Union[str, Path, Dict], structure_only: bool = False, seed: int = 20260923) -> Dict:
    """Checks a contributed specification. `passed` is true when every required check passes."""
    raw = json.loads(Path(source).read_text(encoding="utf-8")) if not isinstance(source, dict) else source
    rows: List[Dict] = []
    try:
        real = fspec.load(raw)
        rows.append(_row("loads", True, "the specification loads", True))
    except fspec.SpecError as err:
        rows.append(_row("loads", False, str(err), True))
        return {"passed": False, "checks": rows}
    prov = raw.get("provenance")
    src = prov.get("source", "") if isinstance(prov, dict) else str(prov or "")
    question = str(raw.get("question", ""))
    field = raw.get("field")
    assumptions = raw.get("assumptions") or []
    rows.append(_row("question", bool(question.strip()) and not _placeholder(question),
                     "states the physical question" if not _placeholder(question) else
                     "replace the template question with the question this material answers", True))
    rows.append(_row("source", bool(src.strip()) and not _placeholder(src),
                     "names where the equations come from" if src.strip() and not _placeholder(src) else
                     "provenance.source must name the paper or book the equations come from", True))
    rows.append(_row("field", isinstance(field, str) and bool(field.strip()) and not _placeholder(field),
                     "names the field of the material" if isinstance(field, str) and not _placeholder(field) else
                     "field must name the field of the material (for example 'soft matter')", True))
    good = [a for a in assumptions if str(a).strip() and not _placeholder(a)]
    rows.append(_row("assumptions", bool(good) and len(good) == len(assumptions),
                     f"{len(good)} assumption(s) stated" if good and len(good) == len(assumptions) else
                     "list the modelling choices a reader could dispute, and remove the template entry", True))
    cited = bool(REFERENCE.search(src)) and not _placeholder(src)
    rows.append(_row("reference", cited,
                     "the source can be looked up (a year, DOI, arXiv identifier or URL)" if cited
                     else "add a citable reference to provenance.source (journal and year, DOI or arXiv)", False))
    if raw.get("kind") == "equations":
        rows.append(_row("control", bool(raw.get("control")),
                         "a control parameter is declared, so write points can be searched along it"
                         if raw.get("control") else
                         "declare the parameter an experiment varies as control, with its range", False))
    for key in ("closure", "observable"):
        value = str(raw.get(key, ""))
        rows.append(_row(key, bool(value.strip()) and not _placeholder(value),
                         f"{key} is stated" if value.strip() and not _placeholder(value) else
                         f"state the {key} (see the glossary of the memory modules)", False))
    from . import predict
    pred = predict.predict(real, np.random.default_rng(seed))
    preds = pred.get("predictions", {})
    summary = {"name": real.name, "structure": {k: v.get("prediction") for k, v in preds.items()
                                                if isinstance(v, dict) and "prediction" in v}}
    if not structure_only:
        from . import discovery
        card = discovery.evaluate(real, np.random.default_rng(seed + 1), quick=True)
        comparison = predict.compare(pred, card)
        summary["card"] = card["verdict"]
        summary["agreement"] = comparison.get("all_consistent")
        summary["mechanism"] = discovery.mechanism_class(card)
        summary["same_mechanism"] = same_mechanism(summary["mechanism"], exclude=real.name)
        rows.append(_row("agreement", bool(comparison.get("all_consistent")),
                         "the calculated card agrees with every structural prediction"
                         if comparison.get("all_consistent") else
                         "a structural prediction disagrees with the calculation: inspect 'comparison' in the card of "
                         "'fieldbridge memory card'", False))
    passed = all(r["passed"] for r in rows if r["required"])
    return {"passed": passed, "checks": rows, "summary": summary}


PREDICTION_NAMES = {"oscillation": "Oscillation", "multistability": "Multistability", "write": "Write point",
                    "holding": "Retention", "lock": "Retention and rewriting times", "network": "Loops"}
VERDICT_NAMES = {"stores": "Stable states", "writes": "Write point", "constructs": "Writing mechanism",
                 "holds": "Retention law", "lock": "Writing protocols"}


def same_mechanism(mechanism: str, exclude: str = "", catalog_json: Union[str, Path] = CATALOG_JSON) -> List[Dict]:
    """Materials of the catalog with the same writing mechanism, from docs/materials.json."""
    path = Path(catalog_json)
    if not path.exists():
        return []
    rows = json.loads(path.read_text(encoding="utf-8")).get("materials", [])
    return [{"name": r["name"], "field": r.get("field", ""), "file": r.get("file", "")}
            for r in rows if r.get("mechanism") == mechanism and r.get("name") != exclude]


def markdown(result: Dict, path: str = "") -> str:
    lines = [f"# Contribution check{': ' + path if path else ''}", "",
             f"**{'Ready to contribute' if result['passed'] else 'Not ready'}**: required checks "
             f"{'all pass' if result['passed'] else 'fail'}.", "", "| Check | Required | Result | Message |",
             "| --- | --- | --- | --- |"]
    for r in result["checks"]:
        lines.append(f"| {r['check']} | {'yes' if r['required'] else 'recommended'} | "
                     f"{'pass' if r['passed'] else 'fail'} | {r['message']} |")
    s = result.get("summary") or {}
    if s.get("structure"):
        lines += ["", "## Predicted from structure", ""] + [f"- **{PREDICTION_NAMES.get(k, k)}**: {v}"
                                                            for k, v in s["structure"].items()]
    if s.get("card"):
        lines += ["", "## Calculated (quick card)", ""] + [f"- **{VERDICT_NAMES.get(k, k)}**: {v}"
                                                          for k, v in s["card"].items()]
    if "same_mechanism" in s:
        others = s["same_mechanism"]
        fields = sorted({o["field"] for o in others if o["field"]})
        lines += ["", "## The same writing mechanism in the catalog", "",
                  f"Mechanism: **{s.get('mechanism')}**. "
                  + (f"{len(others)} material(s) of the catalog, from {len(fields)} field(s), are written the same way:"
                     if others else "No material of the catalog is written this way yet.")]
        lines += [f"- {o['name']} ({o['field']})" for o in others]
    return "\n".join(lines) + "\n"


def specifications(paths: Iterable[Union[str, Path]]) -> List[Path]:
    """Specification files of kind equations or network under the given files and directories (recursively)."""
    out = []
    for p in map(Path, paths):
        files = sorted(p.rglob("*.json")) if p.is_dir() else [p]
        for f in files:
            try:
                kind = json.loads(f.read_text(encoding="utf-8")).get("kind")
            except (json.JSONDecodeError, UnicodeDecodeError):
                continue
            if kind in ("equations", "network"):
                out.append(f)
    return out


def catalog(paths: Iterable[Union[str, Path]], seed: int = 20260923, root: Union[str, Path] = ".") -> Dict:
    """One row per material, with the calculated memory from a quick card."""
    from . import discovery
    rows, skipped = [], []
    for f in specifications(paths):
        raw = json.loads(f.read_text(encoding="utf-8"))
        try:
            real = fspec.load(raw)
        except fspec.SpecError as err:
            skipped.append({"file": str(f), "error": str(err)})
            continue
        card = discovery.evaluate(real, np.random.default_rng(seed), quick=True)
        prov = raw.get("provenance")
        rows.append({"file": str(Path(f).resolve().relative_to(Path(root).resolve())) if root else str(f),
                     "name": real.name, "field": raw.get("field", ""), "question": raw.get("question", ""),
                     "source": prov.get("source", "") if isinstance(prov, dict) else str(prov),
                     "contributor": prov.get("contributor", "") if isinstance(prov, dict) else "",
                     "mechanism": discovery.mechanism_class(card), "verdict": card["verdict"]})
    return {"materials": rows, "skipped": skipped}


def catalog_markdown(result: Dict) -> str:
    def cell(text, n=150):
        text = str(text).replace("|", "\\|").replace("\n", " ")
        return text if len(text) <= n else text[: n - 1].rsplit(" ", 1)[0] + "…"

    rows = sorted(result["materials"], key=lambda r: (str(r["field"]).lower(), r["name"].lower()))
    fields = sorted({str(r["field"]) for r in rows})
    lines = ["# Materials", "",
             f"{len(rows)} materials from {len(fields)} fields. Each row is one specification file. The columns from "
             "the writing mechanism to the retention law are calculated by a quick memory card "
             "(`fieldbridge memory card FILE --quick`). [CONTRIBUTING.md](../CONTRIBUTING.md) describes how a "
             "material is added; this page is regenerated with "
             "`python3 -B -m fieldbridge memory catalog --out docs/materials.md`.", "",
             "| Material | Field | Question | Writing mechanism | Stable states | Write point | Retention law | Source |",
             "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in rows:
        v = r["verdict"]
        link = f"[{cell(r['name'], 60)}](../{r['file']})"
        src = cell(r["source"], 90) + (f" (contributed by {cell(r['contributor'], 40)})" if r.get("contributor") else "")
        lines.append(f"| {link} | {cell(r['field'], 30)} | {cell(r['question'], 110)} | {cell(r['mechanism'], 40)} | "
                     f"{cell(v.get('stores', ''), 50)} | {cell(v.get('writes', ''), 110)} | {cell(v.get('holds', ''), 60)} "
                     f"| {src} |")
    if result.get("skipped"):
        lines += ["", "Not loaded: " + "; ".join(f"`{s['file']}` ({cell(s['error'], 80)})" for s in result["skipped"])]
    lines += ["", "The classifications describe the supplied models. They do not establish that a model describes "
              "the physical material; that needs comparison with experiment and with the literature."]
    return "\n".join(lines) + "\n"
