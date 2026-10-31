"""Discretized Ohmic bath modes instead of Langevin friction: trajectories and the exact linear law.

    python bath_scan.py traj --backend torch --out results/bath_traj.json
    python bath_scan.py law --workers 20 --out results/bath_law.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bath  # noqa: E402
from ptsb import Minimal, law_thermal_speed, linear_law, ratio  # noqa: E402

EPS = 0.04


def points():
    """B: tilt 0.04 at thermal speeds (the ratio hardly depends on friction here: a consistency check only).
    Bf: tilt 0.04 at fast fixed speeds, where a conservative passage departs from the law.
    Bk: curved valley, curvature 0.1 and thermal speeds, where the share depends strongly on friction.
    Bw: as Bk with M = 32 and the band of the modes changed. Each is compared with Langevin trajectories at gamma."""
    P = [dict(id=f"B_g{g}_M{M}_T{kT}", family="B", gamma=g, M=M, kT=kT, n=10000, surface=("minimal",), speed=None,
              w_max=4.0) for g in (0.3, 1.0, 3.0) for M in (1, 2, 4, 8, 16, 32) for kT in (0.01, 0.05)]
    for g in (0.3, 1.0, 3.0):
        for v in (1.0, 3.0):
            for M in (0, 1, 4, 16, 32):         # M = 0: Langevin at gamma, the reference
                P.append(dict(id=f"Bf_g{g}_M{M}_v{v}", family="Bf", gamma=g, M=M, kT=0.005, n=10000,
                              surface=("minimal",), speed=v, w_max=4.0))
        for kT in (0.01, 0.05):
            for M in (0, 1, 2, 4, 8, 16, 32):
                P.append(dict(id=f"Bk_g{g}_M{M}_T{kT}", family="Bk", gamma=g, M=M, kT=kT, n=10000,
                              surface=("curved",), speed=None, w_max=4.0))
            for w_max in (1.0, 16.0):
                P.append(dict(id=f"Bw_g{g}_W{w_max}_T{kT}", family="Bw", gamma=g, M=32, kT=kT, n=10000,
                              surface=("curved",), speed=None, w_max=w_max))
    return P


def build(p):
    from curved import CurvedValley
    return Minimal(EPS) if p["surface"][0] == "minimal" else CurvedValley(kappa=0.1)


def part_traj(args):
    rng = np.random.default_rng(args.seed)
    out = []
    t0 = time.time()
    for p in points():
        s = build(p)
        if p["M"] == 0:                         # the Langevin reference at the same friction, same assignment
            from ptsb import sample_ts
            if p["surface"][0] == "minimal":
                from engine import run_batch
                o = run_batch([dict(surface=s, init=sample_ts(s, p["n"], rng, kT=p["kT"], speed=p["speed"]),
                                    gamma=p["gamma"], kT=p["kT"], t_max=60.0)], backend=args.backend,
                              seed=args.seed)[0][0]
            else:
                from curved import run_batch as run_curved
                o = run_curved([dict(surface=s, init=sample_ts(s, p["n"], rng, kT=p["kT"], speed=p["speed"]),
                                     gamma=p["gamma"], kT=p["kT"], t_max=60.0)], backend=args.backend,
                               seed=args.seed)[0][0]
        else:
            init = bath.sample(s, p["n"], rng, p["kT"], p["gamma"], p["M"], w_max=p["w_max"], speed=p["speed"])
            o = bath.run(s, init, p["gamma"], p["M"], w_max=p["w_max"], backend=args.backend)
        out.append({**p, "outcome": ratio(o)})
        print(p["id"], round(out[-1]["outcome"]["P_plus"], 3), f"{time.time() - t0:.0f} s", flush=True)
        Path(args.out).write_text(json.dumps({"eps": EPS, "seed": args.seed, "points": out}, indent=1))


def law_point(p):
    """The exact law with the modes exists for the straight surface only (family B and Bf)."""
    s = Minimal(EPS)
    t_end = 40.0 if p["gamma"] >= 2 else 12.0
    if p["speed"] is None:
        return p["id"], {"law": bath.law_thermal(s, p["kT"], p["gamma"], p["M"], n_speeds=args_speeds),
                         "markov": law_thermal_speed(s, p["gamma"], p["kT"], t_end=t_end)}
    row = {"markov": {"P_law": linear_law(s, p["speed"], p["gamma"], p["kT"], t_end=t_end)["P_law"]}}
    if p["M"] > 0:
        row["law"] = bath.law(s, p["speed"], p["kT"], p["gamma"], p["M"])
    return p["id"], row


args_speeds = 12


def part_law(args):
    from multiprocessing import Pool
    t0 = time.time()
    P = [p for p in points() if p["family"] in ("B", "Bf")]
    rows = {}
    with Pool(args.workers) as pool:
        for k, row in pool.imap_unordered(law_point, P, chunksize=1):
            rows[k] = row
            print(k, len(rows), f"{time.time() - t0:.0f} s", flush=True)
            Path(args.out).write_text(json.dumps({"laws": rows, "seconds": time.time() - t0}, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("part", choices=["traj", "law"])
    ap.add_argument("--backend", default="numpy")
    ap.add_argument("--seed", type=int, default=20261030)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    (part_traj if args.part == "traj" else part_law)(args)


if __name__ == "__main__":
    main()
