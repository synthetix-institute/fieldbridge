"""Branching on a curved valley, with and without a pseudo-minimum before the VRI: trajectories and laws.

    python curved_scan.py traj --backend torch --only G,G1,H --out results/curved_traj_GH.json   (and I, J)
    python curved_scan.py law --workers 16 --out results/curved_law_all.json

Families (d = 0.25 except J: TS1 at s = 0, the VRI at s = 1/2, TS2 at s = 1; outcome -1 = outside the bend)
  G   curvature alone (eps = 0, symmetric in every static quantity): kappa x friction x temperature, thermal speeds
  G1  curvature alone, gamma = 0, fixed incoming speeds
  H   curvature against a tilt toward the inside of the bend (eps < 0), thermal speeds
  I   a dip of depth a in the profile at s = 0.25 (a shoulder for a < 0.2, a pseudo-minimum with an exit saddle TS1b
      for a >= 0.2), with the curvature kappa = 0.1 or with the tilt eps = 0.02 on a straight valley
  J   an intermediate held between TS1 and an exit saddle TS1b at the energy of TS1 (0, -0.02, -0.04), well 0.1 below
      TS1, the VRI at s = 0.7 after TS1b; same two asymmetries; compared with the same valley without it
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from curved import CurvedValley, curved_law, curved_law_thermal  # noqa: E402
from ptsb import ratio, sample_ts  # noqa: E402

T_BINS = np.r_[0.0, np.geomspace(0.5, 400.0, 25)]


def points():
    P = []
    for kap in (0.02, 0.05, 0.1, 0.2):
        for g in (0.0, 0.3, 1.0, 3.0, 10.0):
            for kT in (0.005, 0.01, 0.02, 0.05, 0.1):
                P.append(dict(id=f"G_k{kap}_g{g}_T{kT}", family="G", surface=dict(kappa=kap), gamma=g, kT=kT,
                              speed=None, n=10000, t_max=80.0 if g >= 3 else 40.0))
    for kap in (0.05, 0.1):
        for kT in (0.005, 0.02):
            for v in (0.3, 1.0, 3.0):
                P.append(dict(id=f"G1_k{kap}_T{kT}_v{v}", family="G1", surface=dict(kappa=kap), gamma=0.0, kT=kT,
                              speed=v, n=10000, t_max=40.0))
    for eps in (-0.01, -0.02, -0.04, -0.08):
        for g in (0.0, 1.0):
            for kT in (0.005, 0.01, 0.02, 0.05, 0.1):
                P.append(dict(id=f"H_e{eps}_g{g}_T{kT}", family="H", surface=dict(kappa=0.1, eps=eps), gamma=g,
                              kT=kT, speed=None, n=10000, t_max=40.0))
    for asym, sf in (("k", dict(kappa=0.1)), ("e", dict(eps=0.02))):
        for a in (0.0, 0.1, 0.15, 0.2, 0.25, 0.3):
            for g in (0.0, 0.3, 1.0, 3.0):
                for kT in (0.02, 0.05):
                    P.append(dict(id=f"I{asym}_a{a}_g{g}_T{kT}", family="I", asym=asym, surface={**sf, "a": a},
                                  gamma=g, kT=kT, speed=None, n=10000, t_max=300.0))
    # an intermediate 0.1 below TS1 between TS1 and an exit saddle TS1b at 0, -0.02 or -0.04 relative to TS1, the VRI
    # at s = 0.7 behind it; "none" is the same valley without the intermediate
    J_PROFILES = {"none": dict(), "b0": dict(a=-0.0447, s_p=0.2, w=0.08, b=0.4028),
                  "b2": dict(a=-0.0493, s_p=0.2, w=0.08, b=0.3805), "b4": dict(a=-0.0549, s_p=0.2, w=0.08, b=0.3576)}
    for asym, sf in (("k", dict(kappa=0.1)), ("e", dict(eps=0.02))):
        for prof, pr in J_PROFILES.items():
            for g in (0.0, 0.1, 0.3, 1.0):
                for kT in (0.02, 0.05):
                    P.append(dict(id=f"J{asym}_{prof}_g{g}_T{kT}", family="J", asym=asym, profile=prof,
                                  surface={**sf, "d": 0.45, **pr}, gamma=g, kT=kT, speed=None, n=10000,
                                  t_max=300.0))
    return P


def part_traj(args):
    from curved import run_batch
    P = points()
    if args.only:
        P = [p for p in P if p["family"] in args.only.split(",")]
    rng = np.random.default_rng(args.seed)
    out = {"dt": args.dt, "seed": args.seed, "backend": args.backend, "t_bins": T_BINS.tolist(), "points": []}
    t0 = time.time()
    # one batch per t_max, so that short runs do not wait for the longest
    for tm in sorted({p["t_max"] for p in P}):
        group = [p for p in P if p["t_max"] == tm]
        batch = []
        for p in group:
            c = CurvedValley(**p["surface"])
            batch.append(dict(surface=c, init=sample_ts(c, p["n"], rng, kT=p["kT"], speed=p["speed"]),
                              gamma=p["gamma"], kT=p["kT"], t_max=tm))
        res = run_batch(batch, backend=args.backend, seed=args.seed + int(tm), dt=args.dt)
        for p, (o, td) in zip(group, res):
            hist = {str(side): np.histogram(td[o == side], T_BINS)[0].tolist() for side in (-1, 1)}
            out["points"].append({**p, "outcome": ratio(o), "time_hist": hist})
        print(f"t_max {tm}: {len(group)} points, {time.time() - t0:.0f} s", flush=True)
    out["seconds"] = time.time() - t0
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1))


def law_point(p):
    c = CurvedValley(**p["surface"])
    t_end = 40.0 if p["gamma"] >= 2 else 12.0
    row = {}
    if p["speed"] is None:
        row["law"] = curved_law_thermal(c, p["gamma"], p["kT"], t_end=t_end)
        row["law_first_order"] = curved_law_thermal(c, p["gamma"], p["kT"], t_end=t_end, closure=False)
        row["law_no_curvature"] = curved_law_thermal(c, p["gamma"], p["kT"], t_end=t_end, curvature=False)
    else:
        for key, cl in (("law", True), ("law_first_order", False)):
            L = curved_law(c, p["speed"], p["gamma"], p["kT"], t_end=t_end, closure=cl)
            row[key] = {k: L[k] for k in ("P_law", "trapped", "v_decision") if k in L}
    if p["family"] in ("I", "J"):
        pp = c.profile_points()
        row["profile"] = {k: (None if v is None else float(v)) for k, v in pp.items()}
        if pp["TS1b"] is not None:
            row["exit_barrier"] = float(c.V0(pp["TS1b"], 0.0) - c.V0(pp["INT"], 0.0))
            row["law_reset"] = curved_law_thermal(c, p["gamma"], p["kT"], t_end=t_end, start="TS1b")
            row["law_reset_first_order"] = curved_law_thermal(c, p["gamma"], p["kT"], t_end=t_end, start="TS1b",
                                                              closure=False)
    return p["id"], row


def part_law(args):
    P = points()
    if args.only:
        P = [p for p in P if p["family"] in args.only.split(",")]
    t0 = time.time()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    rows = {}

    def save():
        Path(args.out).write_text(json.dumps({"seconds": time.time() - t0, "laws": rows}, indent=1))

    if args.workers > 1:
        from multiprocessing import Pool
        with Pool(args.workers) as pool:
            for k, row in pool.imap_unordered(law_point, P, chunksize=1):
                rows[k] = row
                if len(rows) % 10 == 0:
                    save()
                    print(f"{len(rows)} of {len(P)}, {time.time() - t0:.0f} s", flush=True)
    else:
        rows = dict(map(law_point, P))
    save()
    print(f"{len(rows)} laws, {time.time() - t0:.0f} s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("part", choices=["traj", "law"])
    ap.add_argument("--backend", default="numpy")
    ap.add_argument("--only", default="")
    ap.add_argument("--seed", type=int, default=20261030)
    ap.add_argument("--dt", type=float, default=0.005)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    (part_traj if args.part == "traj" else part_law)(args)


if __name__ == "__main__":
    main()
