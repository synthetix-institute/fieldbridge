#!/usr/bin/env python3
"""Replay recorded failure -> verified repair -> reuse through the public CLI, with no external model calls."""
import argparse
import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fieldbridge.cli import main as cli_main
from fieldbridge.memory import learning as ln
from fieldbridge.memory import repair_families as rf


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise SystemExit("Choose a fresh directory")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    task = json.loads((ROOT / "examples/learning/hidden_exchange.json").read_text(encoding="utf-8"))
    registry = args.out_dir / "registry.sqlite"
    families_path = args.out_dir / "families.json"
    evaluations = ln.MoveLibrary(registry)
    ln.retain_record(evaluations, {"spec": task["model"], "targets": {"reduction": {
        "status": "obstructed", "learning_task": task}}}, source="supplied-exchange-equations")
    code = cli_main(["memory", "repair", "--library", str(registry), "--families", str(families_path),
                     "--out-dir", str(args.out_dir / "repair")])
    if code:
        raise SystemExit(code)
    families = rf.FamilyLibrary(families_path, frozen=True)
    digest = families.digest()
    new = copy.deepcopy(task)
    new["model"]["drift"]["x"] = "-a*x+2*b*y"
    for row in new["withheld"]:
        row["initial"] = [v*1.17 for v in row["initial"]]
    warm = rf.construct(new, families, evaluation_library=evaluations)
    cold = ln.learn(new, ln.MoveLibrary(args.out_dir / "cold.sqlite"))
    symbolic = rf.construct(new, strategy="symbolic",
                            evaluation_library=ln.MoveLibrary(args.out_dir / "symbolic.sqlite"))
    assert families.digest() == digest
    assert all(r["status"] == "verified_extension" for r in (warm, cold, symbolic))
    report = {"scope": "supplied deterministic exchange equations; no novelty claim",
              "frozen_library_sha256": digest, "warm": warm, "cold": cold, "symbolic": symbolic}
    (args.out_dir / "campaign_comparison.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k: {f: r[f] for f in ("status", "candidate_evaluations", "seconds")}
                      for k, r in (("warm", warm), ("cold", cold), ("symbolic", symbolic))}, indent=2))


if __name__ == "__main__":
    main()
