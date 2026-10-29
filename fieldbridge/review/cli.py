"""``fieldbridge review``: four checks of the claims of a manuscript.

    fieldbridge review check examples/review/ideal_glass_2404.07492.json --out-dir build/review
    fieldbridge review theorems
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ..core.provenance import provenance, write_report
from . import checks, theorems

BOUNDARY = ("The checks read what the review file states about the manuscript: the systems, the claims, the estimates, "
            "the theorems used and the searches made. A failed check names a specific fault. A passed check is a "
            "necessary condition only; it does not establish that a claim is right, and it does not replace the "
            "calculations and the reading that the review file records.")
HERE = Path(__file__).resolve().parent
LABEL = {"same_system": "Same system", "protocol_independence": "Independence of the protocol",
         "theorem_scope": "Scope of the theorems", "prior_art": "Prior art"}


def render(rep: dict) -> str:
    p = rep["paper"]
    lines = [f"# Review checks: {p.get('title', '')}", ""]
    if p.get("reference"):
        lines += [p["reference"], ""]
    t = rep["tally"]
    lines += [f"{len(rep['claims'])} claims; checks passed {t['pass']}, warned {t['warn']}, failed {t['fail']}, "
              f"not applicable {t[checks.NA]}.", ""]
    for c in rep["claims"]:
        lines += [f"## {c['id']}: {c['statement']}", ""]
        for r in c["results"]:
            lines += [f"- **{LABEL[r['check']]}: {r['status']}.** {r['explanation']}"]
        lines += [""]
    lines += ["## Scope", "", BOUNDARY, ""]
    return "\n".join(lines)


def cmd_check(args) -> int:
    doc = json.loads(Path(args.review).read_text(encoding="utf-8"))
    if doc.get("schema") != "fieldbridge-review/1":
        raise SystemExit("not a fieldbridge-review/1 file")
    rep = checks.review(doc)
    rep["provenance"] = provenance(HERE, BOUNDARY, doc, libraries=("sympy",))
    md = render(rep)
    write_report(Path(args.out_dir), "review", rep, md)
    print(md)
    return 1 if args.strict and rep["tally"]["fail"] else 0


def cmd_theorems(args) -> int:
    for name, statement in theorems.table():
        print(f"{name}: {statement}")
    return 0


def add_parser(sub) -> None:
    rv = sub.add_parser("review", help="Four checks of the claims of a manuscript: same system, independence of the "
                                       "protocol, scope of the theorems, prior art.")
    s = rv.add_subparsers(required=True)
    c = s.add_parser("check", help="Run the four checks on a review file (schema fieldbridge-review/1).")
    c.add_argument("review")
    c.add_argument("--out-dir", default="build/review")
    c.add_argument("--strict", action="store_true", help="exit with status 1 when a check fails")
    c.set_defaults(func=cmd_check)
    t = s.add_parser("theorems", help="List the theorems whose hypotheses the scope check knows.")
    t.set_defaults(func=cmd_theorems)
