"""Mechanisms of a quantum-computing hardware model: the three derivations of the memory language with their laws,
the law at every obstruction (memory complete), and the single changes of the model (memory changes).

Each model is the classical image of an element of quantum-computing hardware: the mean field of its Lindbladian with
the vacuum noise of its loss channels, with published parameters (the source is in its provenance). One process per
model; the results are written to OUT/<model>.json and collected by report.py.

    python3 -B examples/hardware/survey.py examples/hardware/kerr_cat_kpo.json --out examples/hardware/results

The full survey of the 31 models takes about two hours on 16 cores; the laws at the obstructions take most of it.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fieldbridge.memory import changes, codiscovery, completion, spec  # noqa: E402

SEED = 20260923
TARGETS = ("symmetric-write", "threshold-write", "phase-locking")


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items() if k not in ("final", "spec")}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


def survey(path: Path, trajectories: int = 400, replicates: int = 8, samples: int = 8000, certificate: int = 4000,
           budget: int = 60, save=None) -> dict:
    """The three derivations, the laws at the obstructions and the single changes of one model; `save` is called
    with the partial result after each stage."""
    real = spec.load(path)
    out = {"model": path.stem, "name": real.spec.get("name"), "provenance": real.spec.get("provenance"),
           "derivations": {}, "seconds": {}}
    for t in TARGETS:
        t0 = time.time()
        rep = codiscovery.codiscover([real], SEED, target=t, check_law=True, n_traj=trajectories, replicates=replicates)
        row = rep["rows"][0]
        out["derivations"][t] = {k: row.get(k) for k in ("status", "obstruction", "obstruction_short", "class", "word",
                                                         "law_constant")}
        out["seconds"][t] = round(time.time() - t0, 1)
        if save:
            save(clean(out))
    t0 = time.time()
    out["completion"] = completion.complete(real, np.random.default_rng(1), n_samples=samples, n_traj=certificate,
                                            floors=(0.5, 1.0, 2.0))
    out["seconds"]["completion"] = round(time.time() - t0, 1)
    if save:
        save(clean(out))
    t0 = time.time()
    out["changes"] = changes.changes(real.spec, budget=budget, seed=0)
    out["seconds"]["changes"] = round(time.time() - t0, 1)
    return clean(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("spec", type=Path)
    ap.add_argument("--out", type=Path, default=ROOT / "examples/hardware/results")
    ap.add_argument("--trajectories", type=int, default=400, help="Trajectories per replicate of a law constant.")
    ap.add_argument("--replicates", type=int, default=8)
    ap.add_argument("--samples", type=int, default=8000, help="Samples of the linear stage of a law at an obstruction.")
    ap.add_argument("--certificate", type=int, default=4000, help="Trajectories of its certificate per bias.")
    ap.add_argument("--budget", type=int, default=60, help="Time limit of a derivation after a change, in seconds.")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    target = args.out / f"{args.spec.stem}.json"
    partial = args.out / f"{args.spec.stem}.partial.json"

    def save(r):
        partial.write_text(json.dumps(r, indent=1, default=float) + "\n", encoding="utf-8")

    res = survey(args.spec, args.trajectories, args.replicates, args.samples, args.certificate, args.budget, save)
    target.write_text(json.dumps(res, indent=1, default=float) + "\n", encoding="utf-8")
    if partial.exists():
        partial.unlink()
    print(args.spec.stem, "done", res["seconds"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
