"""Friction from cubic coupling to further modes: trajectories, the equivalent harmonic system and the law.

    python anharm_scan.py traj --backend torch --part 0/4 --out results/anharm_traj_0.json     (parts 1/4 to 3/4)
    python anharm_scan.py merge --out results/anharm_traj.json     (joins results/anharm_traj_?.json, keeps what is used)
    python anharm_scan.py law --workers 16 --out results/anharm_law.json

Surfaces, chosen because the share depends on friction there (the first two), or as a consistency check (the third):
  curved  curved valley, curvature 0.1, thermal speeds, kT = 0.01 and 0.05
  fast    tilt 0.04, speed 3, kT = 0.005
  tilt    tilt 0.04, thermal speeds, kT = 0.01 and 0.05
Families (three random tables each, M = 16 modes with frequencies 0.5 to 8)
  F   each cubic form alone, g = 0.1, 0.3, 1, 3 (quadratic forms 0.1, 0.3, 0.5, 1, where a table stays bound); the sets
      disp (xxQ, xyQ, yyQ), quad (xQQ, yQQ), all, all + QQQ; in the mixed sets the quadratic part stops at g = 0.3
  E   the equivalent harmonic system of every point of F with a quadratic form
  G   the limit in which the displaced modes follow exactly (inertia, no friction, no noise), every displacement point of F
  M   the set all at g = 1 with 4, 8 and 32 modes
  Q   the set all at g = 1 and kT = 0.03 with classical modes, with the Bose populations and with zero-point energy
  A   no tilt, modes stiffening on one side (B_jj = g w_j), g = 0.1 and 0.3, kT = 0.005 to 0.1
  R   references: the surfaces with no coupling, and harmonic modes with a counterterm at known friction
  S   the forms xxQ, xyQ, disp at g = 1 on the curved valley with the frequencies multiplied by 2, 4, 8 at fixed displacement
      (which must converge to G), and with only the modes below or above frequency 2
  D   harmonic modes of known friction on x alone or on y alone, on the curved valley and on the straight tilted valley; on
      the curved valley also on s alone or on n alone, along the path and across it
  P   the forms xxQ, xyQ, disp on the curved valley with the coupling written in (s, n) instead of (x, y)
  Z   the two limits of S by band: the modes below frequency 2 frozen, and the modes above frequency 2 following exactly
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import zlib
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import anharm as an  # noqa: E402
import bath  # noqa: E402
from curved import CurvedValley  # noqa: E402
from ptsb import Minimal, ratio, sample_ts  # noqa: E402

SETS = {"xxQ": ("xxQ",), "xyQ": ("xyQ",), "yyQ": ("yyQ",), "xQQ": ("xQQ",), "yQQ": ("yQQ",),
        "disp": ("xxQ", "xyQ", "yyQ"), "quad": ("xQQ", "yQQ"), "all": an.FORMS, "allQ3": an.FORMS + ("QQQ",)}
QUAD = ("xQQ", "yQQ", "quad")
GS = (0.1, 0.3, 1.0, 3.0)
GS_QUAD = (0.1, 0.3, 0.5, 1.0)                    # the quadratic tables lose binding beyond about 0.5 (16 modes)
SEEDS = (1, 2, 3)
GQ_MAX = 0.3                                      # quadratic tables in the mixed sets (they lose binding beyond it)
CONDITIONS = [("curved", 0.01, None), ("curved", 0.05, None), ("fast", 0.005, 3.0), ("tilt", 0.01, None),
              ("tilt", 0.05, None)]


def surface(name):
    return {"curved": lambda: CurvedValley(kappa=0.1), "fast": lambda: Minimal(0.04), "tilt": lambda: Minimal(0.04),
            "flat": lambda: Minimal(0.0)}[name]()


def table(p):
    """The coupling table of a point (None if no table keeps every mode bound at that strength)."""
    if p["family"] == "D":
        w, c = bath.modes(p["gamma"], p["M"], 4.0)
        K0 = float((c ** 2 / w ** 2).sum())
        frame = "path" if p["axis"] in "sn" else "xy"
        return (an.Tables(w=w, Lx=c, Kx=K0, frame=frame) if p["axis"] in "xs" else
                an.Tables(w=w, Ly=c, Ky=K0, frame=frame)), 0
    if p["family"] == "P":
        T, r = an.draw(SETS[p["set"]], p["g"], p["M"], p["seed"])
        T.frame = "path"
        return T, r
    if p["family"] in ("S", "Z", "ZG"):
        T, r = an.draw(SETS[p["set"]], p["g"], p["M"], p["seed"])
        if p.get("band"):
            keep = T.w <= 2.0 if p["band"] == "low" else T.w > 2.0
            T.dxx, T.dxy, T.dyy = T.dxx * keep, T.dxy * keep, T.dyy * keep
        T.w = T.w * p.get("lam", 1.0)
        return T, r
    if p["family"] == "A":
        return an.entropic(p["g"], p["M"]), 0
    if p["family"] == "R":
        if p["gamma"] == 0:
            return an.Tables(w=an.frequencies(p["M"])), 0
        w, c = bath.modes(p["gamma"], p["M"], 4.0)
        K0 = float((c ** 2 / w ** 2).sum())
        return an.Tables(w=np.r_[w, w], Lx=np.r_[c, 0 * c], Ly=np.r_[0 * c, c], Kx=K0, Ky=K0), 0
    forms = SETS[p["set"]]
    if p["set"] in ("disp", "all", "allQ3") and p["g"] > GQ_MAX:
        Td, r = an.draw(tuple(f for f in forms if f not in ("xQQ", "yQQ", "QQQ")), p["g"], p["M"], p["seed"])
        rest = tuple(f for f in forms if f in ("xQQ", "yQQ", "QQQ"))
        if not rest:
            return Td, r
        Tq, r2 = an.draw(rest, GQ_MAX, p["M"], p["seed"])
        if Tq is None:
            return None, r2
        Td.A, Td.B, Td.tri = Tq.A, Tq.B, Tq.tri
        return Td, r + r2
    return an.draw(forms, p["g"], p["M"], p["seed"])


def points():
    P = []
    for name, kT, v in CONDITIONS:
        tag = f"{name}_T{kT}"
        P.append(dict(id=f"R_{tag}", family="R", name=name, kT=kT, speed=v, gamma=0.0, M=16, n=20000))
        for g in (0.3, 1.0, 3.0):
            P.append(dict(id=f"R_{tag}_g{g}", family="R", name=name, kT=kT, speed=v, gamma=g, M=32, n=10000))
        for st in SETS:
            for g in (GS_QUAD if st in QUAD else GS):
                for seed in SEEDS:
                    base = dict(name=name, kT=kT, speed=v, set=st, g=g, seed=seed, M=16, n=10000)
                    P.append(dict(id=f"F_{tag}_{st}_g{g}_s{seed}", family="F", **base))
                    if st in QUAD:
                        P.append(dict(id=f"E_{tag}_{st}_g{g}_s{seed}", family="E", **base))
                    if st in ("xxQ", "xyQ", "yyQ", "disp"):
                        P.append(dict(id=f"G_{tag}_{st}_g{g}_s{seed}", family="G", **base))
        for M in (4, 8, 32):
            for seed in SEEDS:
                P.append(dict(id=f"M_{tag}_M{M}_s{seed}", family="M", name=name, kT=kT, speed=v, set="all", g=1.0,
                              seed=seed, M=M, n=10000))
    for name in ("curved", "tilt"):
        for occ in ("classical", "bose", "wigner"):
            for seed in SEEDS:
                P.append(dict(id=f"Q_{name}_{occ}_s{seed}", family="Q", name=name, kT=0.03, speed=None, set="all",
                              g=1.0, seed=seed, M=16, n=10000, occupation=occ))
    for g in (0.1, 0.3):
        for kT in (0.005, 0.01, 0.02, 0.05, 0.1):
            P.append(dict(id=f"A_g{g}_T{kT}", family="A", name="flat", kT=kT, speed=None, g=g, M=16, n=20000))
    for kT in (0.01, 0.05):
        for st in ("xxQ", "xyQ", "disp"):
            for seed in SEEDS:
                base = dict(family="S", name="curved", kT=kT, speed=None, set=st, g=1.0, seed=seed, M=16, n=10000)
                for lam in (2.0, 4.0, 8.0):
                    P.append(dict(id=f"S_T{kT}_{st}_lam{lam}_s{seed}", lam=lam, dt=0.00125 if lam > 2 else 0.0025, **base))
                for band in ("low", "high"):
                    P.append(dict(id=f"S_T{kT}_{st}_{band}_s{seed}", band=band, **base))
        for axis in ("x", "y"):
            for g in (0.3, 1.0, 3.0):
                P.append(dict(id=f"D_T{kT}_{axis}_g{g}", family="D", name="curved", kT=kT, speed=None, axis=axis, gamma=g,
                              M=32, n=10000))
    for kT in (0.01, 0.05):
        for g in (0.3, 1.0, 3.0):
            for axis in ("x", "y"):
                P.append(dict(id=f"Dt_T{kT}_{axis}_g{g}", family="D", name="tilt", kT=kT, speed=None, axis=axis, gamma=g,
                              M=32, n=10000))
            for axis in ("s", "n"):
                P.append(dict(id=f"Dp_T{kT}_{axis}_g{g}", family="D", name="curved", kT=kT, speed=None, axis=axis, gamma=g,
                              M=32, n=10000))
        for st in ("xxQ", "xyQ", "disp"):
            for g in (0.3, 1.0, 3.0):
                for seed in SEEDS:
                    P.append(dict(id=f"P_T{kT}_{st}_g{g}_s{seed}", family="P", name="curved", kT=kT, speed=None, set=st, g=g,
                                  seed=seed, M=16, n=10000))
            for seed in SEEDS:
                base = dict(name="curved", kT=kT, speed=None, set=st, g=1.0, seed=seed, M=16, n=10000)
                P.append(dict(id=f"Z_T{kT}_{st}_s{seed}", family="Z", band="low", **base))
                P.append(dict(id=f"ZG_T{kT}_{st}_s{seed}", family="ZG", band="high", **base))
    return P


def part_traj(args):
    k, m = (int(v) for v in args.part.split("/"))
    P = [p for i, p in enumerate(points()) if i % m == k]
    if args.only:
        P = [p for p in P if p["family"] in args.only.split(",")]
    if args.ids:
        P = [p for p in P if p["id"].startswith(tuple(args.ids.split(",")))]
    if args.limit:
        P = P[:args.limit]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out = []
    t0 = time.time()
    for p in P:
        T, redraws = table(p)
        row = dict(p, redraws=redraws)
        if T is None:
            row["unbound"] = True
            out.append(row)
            continue
        row["scale"] = an.friction_scale(T, p["kT"])
        if p["family"] == "E":
            T = an.equivalent(T, p["kT"])
        s = surface(p["name"])
        rng = np.random.default_rng([args.seed, zlib.crc32(p["id"].encode())])
        if p["family"] in ("G", "ZG"):
            r = an.run_following(s, T, sample_ts(s, p["n"], rng, kT=p["kT"], speed=p["speed"]), backend=args.backend)
            row["outcome"] = ratio(r["outcome"])
            row["energy_drift_rms"] = r["energy_drift_rms"]
            row["t_mean"] = float(np.nanmean(r["time"]))
            out.append(row)
            print(p["id"], round(row["outcome"]["P_plus"], 3), f"{time.time() - t0:.0f} s", flush=True)
            continue
        init = an.sample(s, T, p["n"], rng, p["kT"], speed=p["speed"], occupation=p.get("occupation", "classical"),
                         backend=args.backend)
        r = an.run(s, T, init, dt=p.get("dt", 0.005), backend=args.backend, frozen=p["family"] == "Z")
        row["outcome"] = ratio(r["outcome"])
        row["bad"] = int(np.sum(r["outcome"] == 3))
        row["accept"] = init["accept"]
        row["gamma_energy"] = an.gamma_from_energy(r, p["kT"]) if "int_p2" in r else None
        row.update({k_: r[k_] for k_ in ("dE_surface", "dE_coupling", "dE_modes", "dE_left", "dE_left_by_mode",
                                         "E_following", "int_p2", "t_mean", "energy_drift_rms") if k_ in r})
        row["w"] = T.w.tolist()
        out.append(row)
        print(p["id"], round(row["outcome"]["P_plus"], 3), row["gamma_energy"], f"{time.time() - t0:.0f} s", flush=True)
        if len(out) % 5 == 0:
            Path(args.out).write_text(json.dumps({"seed": args.seed, "points": out}, indent=0))
    Path(args.out).write_text(json.dumps({"seed": args.seed, "points": out, "seconds": time.time() - t0}, indent=0))


def law_point(p):
    T, _ = table(p)
    if T is None:
        return p["id"], {"unbound": True}
    s = surface(p["name"])
    if p["speed"] is None:
        return p["id"], an.law_thermal(s, T, p["kT"])
    L = an.law(s, T, p["speed"], p["kT"])
    return p["id"], {k: L[k] for k in ("P_law", "P_law_mean_force", "x_decision") if k in L}


def part_law(args):
    P = [p for p in points() if p["family"] in ("F", "A", "M") and p["name"] != "curved"
         and (p["family"] != "F" or p["set"] != "allQ3")]
    if args.limit:
        P = P[:args.limit]
    t0 = time.time()
    rows = {}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    def save():
        Path(args.out).write_text(json.dumps({"seconds": time.time() - t0, "laws": rows}, indent=0))

    from multiprocessing import Pool
    with Pool(args.workers) as pool:
        for k, row in pool.imap_unordered(law_point, P, chunksize=1):
            rows[k] = row
            if len(rows) % 10 == 0:
                save()
                print(f"{len(rows)} of {len(P)}, {time.time() - t0:.0f} s", flush=True)
    save()
    print(f"{len(rows)} laws, {time.time() - t0:.0f} s", flush=True)


def part_merge(args):
    out = Path(args.out)
    files = sorted(out.parent.glob(out.stem + "_?.json"))
    drop = ("w", "dE_left_by_mode", "gamma_energy", "dE_left", "dE_modes", "dE_coupling")
    P = [{k: v for k, v in p.items() if k not in drop} for f in files for p in json.loads(f.read_text())["points"]]
    P.sort(key=lambda p: p["id"])
    out.write_text(json.dumps({"seed": args.seed, "points": P}, indent=0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("part_name", choices=["traj", "law", "merge"])
    ap.add_argument("--backend", default="numpy")
    ap.add_argument("--part", default="0/1")
    ap.add_argument("--only", default="")
    ap.add_argument("--ids", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--seed", type=int, default=20261030)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    {"traj": part_traj, "law": part_law, "merge": part_merge}[args.part_name](args)


if __name__ == "__main__":
    main()
