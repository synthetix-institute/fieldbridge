#!/usr/bin/env python3
"""Small reproducible learning experiment; no corpus sweep or external model call."""
import argparse
import copy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fieldbridge.memory.learning import MoveLibrary, learn


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    library_path = out / "moves.sqlite"
    if library_path.exists():
        raise SystemExit("Choose a fresh directory: the cold/warm comparison needs an empty library")
    library = MoveLibrary(library_path)
    examples = ROOT / "examples/learning"
    first = json.loads((examples / "quintic_saturation.json").read_text(encoding="utf-8"))
    transfer = copy.deepcopy(first)
    transfer["model"]["parameters"]["c"] = .2
    transfer["model"]["drift"]["x"] = "a*x-c*x**3-2*b*x**5"
    transfer["candidate_drift"]["u"] = "a*u-c*u**3"
    records = {"quintic_cold": learn(first, library), "transfer_warm": learn(transfer, library),
               "transfer_without_library": learn(transfer)}
    hidden = json.loads((examples / "hidden_exchange.json").read_text(encoding="utf-8"))
    records["hidden_initial_condition"] = learn(hidden, library)
    square = copy.deepcopy(hidden)
    square["model"]["drift"]["x"] = "-a*x+b*y**2"
    records["derived_squared_observable"] = learn(square, library)
    records["moving_branch"] = learn(json.loads((examples / "moving_branch.json").read_text(encoding="utf-8")), library)
    summary = {k: {"status": r["status"], "candidate_evaluations": r["candidate_evaluations"],
                   "selection_origin": r.get("selection_origin"), "coordinates": r.get("candidate", {}).get("coordinates"),
                   "drift": r.get("candidate", {}).get("drift"),
                   "withheld_error": r.get("withheld", {}).get("max_error"),
                   "omission_error": r.get("withheld", {}).get("max_omission_error")}
               for k, r in records.items()}
    report = {"scope": "Authored equation-level verification cases; no claim of new physics or source recovery",
              "runs": records, "summary": summary}
    (out / "experiment.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if any(r["status"] != "verified_extension" for r in records.values()):
        raise SystemExit(1)
    if not records["transfer_warm"]["candidate_evaluations"] < records["transfer_without_library"]["candidate_evaluations"]:
        raise SystemExit("No improvement over the empty-library search")


if __name__ == "__main__":
    main()
