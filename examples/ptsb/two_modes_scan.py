"""Two transverse modes along one path: four products (or three), trajectories and the law.

    python two_modes_scan.py traj --backend torch --out results/two_modes_traj.json
    python two_modes_scan.py law --workers 16 --out results/two_modes_law.json
    python two_modes_scan.py sep --workers 16 --out results/two_modes_sep.json

Families
  T  c = 0.5, products in the four quadrants: signs of the two modes against the product of two sign laws
  W  c = 1.5, products on the axes: which mode wins against P(|y| > |z|), over the delay D of the second instability
     and a tilt eps_z that removes one z product (three products)
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from two_modes import TwoModes, law_thermal, run, sample, separatrix_thermal, shares  # noqa: E402


def points():
    P = []
    for g in (0.0, 1.0, 3.0):
        for kT in (0.01, 0.05):
            for ey, ez in ((0.02, 0.0), (0.02, 0.04)):
                for D in (0.0, 0.1, 0.3):
                    P.append(dict(id=f"T_g{g}_T{kT}_ez{ez}_D{D}", family="T", gamma=g, kT=kT,
                                  surface=dict(c=0.5, D=D, eps_y=ey, eps_z=ez)))
            for ez in (0.0, 0.05, 0.2):
                for D in (0.0, 0.02, 0.05, 0.1, 0.2):
                    P.append(dict(id=f"W_g{g}_T{kT}_ez{ez}_D{D}", family="W", gamma=g, kT=kT,
                                  surface=dict(c=1.5, D=D, eps_y=0.0, eps_z=ez)))
    return P


def part_traj(args):
    rng = np.random.default_rng(args.seed)
    rows, t0 = [], time.time()
    for i, p in enumerate(points()):
        s = TwoModes(**p["surface"])
        o = run(s, sample(s, args.n, rng, p["kT"]), p["gamma"], p["kT"], t_max=80.0 if p["gamma"] >= 3 else 40.0,
                backend=args.backend, seed=args.seed + i)
        rows.append({**p, "n": args.n, "result": shares(s, *o)})
        print(p["id"], f"{time.time() - t0:.0f} s", flush=True)
        Path(args.out).write_text(json.dumps({"seed": args.seed, "points": rows}, indent=1))


def law_point(p):
    return p["id"], law_thermal(TwoModes(**p["surface"]), p["gamma"], p["kT"], t_end=40.0 if p["gamma"] >= 2 else 14.0)


def part_law(args):
    from multiprocessing import Pool
    rows, t0 = {}, time.time()
    with Pool(args.workers) as pool:
        for k, row in pool.imap_unordered(law_point, points(), chunksize=1):
            rows[k] = row
            if len(rows) % 10 == 0:
                print(len(rows), f"{time.time() - t0:.0f} s", flush=True)
                Path(args.out).write_text(json.dumps({"laws": rows}, indent=1))
    Path(args.out).write_text(json.dumps({"laws": rows, "seconds": time.time() - t0}, indent=1))
    print(len(rows), "laws,", f"{time.time() - t0:.0f} s", flush=True)


def sep_point(p):
    return p["id"], separatrix_thermal(TwoModes(**p["surface"]), p["gamma"], p["kT"],
                                       t_end=40.0 if p["gamma"] >= 2 else 14.0)


def part_sep(args):
    from multiprocessing import Pool
    P = [p for p in points() if p["family"] == "W"]
    rows, t0 = {}, time.time()
    with Pool(args.workers) as pool:
        for k, row in pool.imap_unordered(sep_point, P, chunksize=1):
            rows[k] = row
            print(len(rows), f"{time.time() - t0:.0f} s", flush=True)
            Path(args.out).write_text(json.dumps({"laws": rows}, indent=1))
    Path(args.out).write_text(json.dumps({"laws": rows, "seconds": time.time() - t0}, indent=1))
    print(len(rows), "laws,", f"{time.time() - t0:.0f} s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("part", choices=["traj", "law", "sep"])
    ap.add_argument("--backend", default="numpy")
    ap.add_argument("--seed", type=int, default=20261030)
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    {"traj": part_traj, "law": part_law, "sep": part_sep}[args.part](args)


if __name__ == "__main__":
    main()
