"""Products on the axes: what becomes of a trajectory after its first entry into a product region.

    python two_modes_fate.py --out results/two_modes_fate.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from two_modes import TwoModes, fate  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/two_modes_fate.json")
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    rows, t0 = [], time.time()
    for D in (0.05, 0.1, 0.2):
        for gamma in (1.0, 3.0):
            for kT in (0.01, 0.05):
                times = (5, 10, 20, 40, 60, 120) if gamma >= 3 else (5, 10, 20, 40, 60)
                r = fate(TwoModes(c=1.5, D=D), gamma, kT, n=args.n, times=times, seed=args.seed)
                rows.append({"D": D, "gamma": gamma, "kT": kT, **r})
                print(f"D {D} gamma {gamma} kT {kT}: first entry {r['z_share_first_entry']:.3f}, later {r['z_share']}, "
                      f"well {r['z_well_over_kT']:.1f} kT ({time.time() - t0:.0f} s)", flush=True)
                Path(args.out).write_text(json.dumps({"seed": args.seed, "points": rows}, indent=1))


if __name__ == "__main__":
    main()
