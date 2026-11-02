#!/usr/bin/env python3
"""Bounded regression, recorded-repair and frozen-transfer validation; writes only to a fresh output folder."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.out_dir.resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit("Choose a fresh output directory")
    out.mkdir(parents=True, exist_ok=True)
    files = [*ROOT.glob("fieldbridge/memory/*.py"), *ROOT.glob("tests/test_*learning.py"),
             ROOT / "tests/test_repair_families.py", ROOT / "tests/test_periodic_construction.py",
             *ROOT.glob("scripts/check_repair*.py")]
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    manifest = {"schema": "fieldbridge-constructor-validation/1", "python": sys.version,
                "implementation": hashes, "stages": [], "completed": False}
    commands = [
        ("regressions", ["-m", "pytest", "-q", "-p", "no:cacheprovider", "tests/test_constructor_learning.py",
                         "tests/test_repair_families.py", "tests/test_periodic_construction.py"]),
        ("campaign", ["scripts/check_repair_campaign.py", "--out-dir", str(out / "campaign")]),
        ("transfer", ["scripts/check_repair_transfer.py", "--out-dir", str(out / "transfer")]),
    ]
    env = {**os.environ, "PYTHONPATH": str(ROOT), "PYTHONDONTWRITEBYTECODE": "1",
           "OMP_NUM_THREADS": "2", "OPENBLAS_NUM_THREADS": "2", "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1"}
    for name, command in commands:
        started = time.monotonic()
        with (out / (name + ".log")).open("w") as log:
            process = subprocess.run([sys.executable, "-B", *command], cwd=ROOT, env=env,
                                     stdout=log, stderr=subprocess.STDOUT)
        manifest["stages"].append({"name": name, "command": command, "exit_code": process.returncode,
                                   "seconds": time.monotonic()-started})
        manifest["completed"] = len(manifest["stages"]) == len(commands) and process.returncode == 0
        (out / "validation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        print(json.dumps(manifest["stages"][-1]), flush=True)
        if process.returncode:
            raise SystemExit(process.returncode)
    if any(hashlib.sha256((ROOT / p).read_bytes()).hexdigest() != value for p, value in hashes.items()):
        manifest["completed"] = False
        manifest["error"] = "Implementation changed during validation"
        (out / "validation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        raise SystemExit("Implementation changed during validation")


if __name__ == "__main__":
    main()
