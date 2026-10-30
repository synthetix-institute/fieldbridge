"""Branching ratios on the model surfaces: trajectories and the linear-stage law over asymmetry, temperature, friction,
incoming energy and the sampling of the zero-point motion.

    python scan.py traj --backend torch --out results/traj.json     # trajectories (GPU if available)
    python scan.py law --out results/law.json                       # the law from the surfaces (no trajectories)

Both parts build the same list of points (`points()`); each point has an id.

Families
  A  Minimal surface, tilt eps: friction gamma x temperature kT x eps, thermal (flux-weighted) speeds
  A1 Minimal surface, tilt eps, gamma = 0, fixed incoming speeds
  B  Minimal surface, cubic asymmetry c (no force on the line y = 0): the law predicts 1/2
  C  Minimal surface, tilt eps, zero-point motion only (kT = 0): quasiclassical normal-mode sampling (fixed energy,
     uniform phase) against the Gaussian Wigner distribution of the same mean energy, without friction and with
     friction and no noise
  D  Collins et al. (2013) surface, xstar in {2.00, 2.05, 2.10}: thermal sampling over gamma and kT
  E  Collins et al. (2013) protocol: fixed energy E above TS1 on the dividing line, damping gamma without noise,
     product assigned by the well occupied at the end (their definition) and at first entry
  F  Minimal surface, tilt eps = 0.02, the VRI moved through the TS line (d from -0.45 to 0.15; for d < -1/4 TS1
     splits into one saddle per product): canonical flux ensemble on the line x = 0 with the full Boltzmann weight;
     compared with the law (VRI after the TS line) and with the transition-state split (VRI before it)
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ptsb import (Collins, Minimal, law_dividing_line, law_thermal_speed, linear_law, qct_probability,  # noqa: E402
                  ratio, sample_dividing_line, sample_ts, ts_frame, tst_split)

HBAR = 0.04          # zero-point energy hbar omega / 2 = 0.02 at the TS1 transverse frequency omega = 1


def points():
    P = []
    for g in (0.0, 0.3, 1.0, 3.0, 10.0):
        for kT in (0.005, 0.01, 0.02, 0.05, 0.1):
            for eps in (0.005, 0.01, 0.02, 0.04, 0.08):
                P.append(dict(id=f"A_g{g}_T{kT}_e{eps}", family="A", surface=("minimal", eps, 0.0), gamma=g, kT=kT,
                              speed=None, transverse="thermal", n=10000, t_max=80.0 if g >= 3 else 40.0))
    for kT in (0.005, 0.02, 0.1):
        for eps in (0.01, 0.04):
            for v in (0.3, 1.0, 3.0):
                P.append(dict(id=f"A1_T{kT}_e{eps}_v{v}", family="A1", surface=("minimal", eps, 0.0), gamma=0.0,
                              kT=kT, speed=v, transverse="thermal", n=10000, t_max=40.0))
    for g in (0.0, 3.0):
        for kT in (0.005, 0.02, 0.05, 0.1):
            for c in (0.1, 0.2, 0.4):
                P.append(dict(id=f"B_g{g}_T{kT}_c{c}", family="B", surface=("minimal", 0.0, c), gamma=g, kT=kT,
                              speed=None, transverse="thermal", n=10000, t_max=80.0 if g >= 3 else 40.0))
    for g in (0.0, 3.0):
        for tr in ("qct", "wigner"):
            for eps in (0.0025, 0.005, 0.01, 0.02, 0.04):
                P.append(dict(id=f"C_g{g}_{tr}_e{eps}", family="C", surface=("minimal", eps, 0.0), gamma=g, kT=0.0,
                              speed=1.0, transverse=tr, hbar=HBAR, n=10000, t_max=80.0 if g >= 3 else 40.0))
    for xs in (2.0, 2.05, 2.1):
        for g in (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0):
            for kT in (0.01, 0.05, 0.2):
                P.append(dict(id=f"D_x{xs}_g{g}_T{kT}", family="D", surface=("collins", xs), gamma=g, kT=kT,
                              speed=None, transverse="thermal", n=10000, t_max=80.0 if g >= 2 else 40.0))
    for xs in (2.0, 2.05, 2.1):
        for E in (0.01, 0.1):
            for g in (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 1.0):
                P.append(dict(id=f"E_x{xs}_E{E}_g{g}", family="E", surface=("collins", xs), gamma=g, kT=0.0,
                              energy=E, n=5000, t_max=120.0, rule="final"))
    for d in (-0.45, -0.4, -0.35, -0.3, -0.27, -0.25, -0.22, -0.15, -0.05, 0.05, 0.15):
        for kT in (0.01, 0.03, 0.1):
            for g in (0.0, 3.0):
                P.append(dict(id=f"F_d{d}_T{kT}_g{g}", family="F", surface=("minimal", 0.02, 0.0, d), gamma=g, kT=kT,
                              sampling="dividing line", n=10000, t_max=80.0 if g >= 3 else 40.0))
    return P


def build(spec):
    return Minimal(*spec[1:]) if spec[0] == "minimal" else Collins(spec[1])


def sample_microcanonical(s, n, rng, E):
    """Fixed total energy E above TS1 on the dividing line, flux weighted: the transverse energy is uniform on
    [0, E] (one transverse mode, harmonic at TS1), its phase uniform, the rest is kinetic energy along u."""
    q0, u, w, om, _ = ts_frame(s)
    Ew = rng.uniform(0.0, E, n)
    ph = rng.uniform(0.0, 2 * np.pi, n)
    qw, pw = np.sqrt(2 * Ew) / om * np.cos(ph), -np.sqrt(2 * Ew) * np.sin(ph)
    pu = np.sqrt(2 * (E - Ew))
    return (q0[0] + qw * w[0], q0[1] + qw * w[1], pu * u[0] + pw * w[0], pu * u[1] + pw * w[1])


def part_traj(args):
    from engine import run_batch
    P = points()
    if args.only:
        P = [p for p in P if p["family"] in args.only.split(",")]
    rng = np.random.default_rng(args.seed)
    batch = []
    for p in P:
        s = build(p["surface"])
        if p["family"] == "E":
            init = sample_microcanonical(s, p["n"], rng, p["energy"])
        elif p["family"] == "F":
            init = sample_dividing_line(s, p["n"], rng, p["kT"])
        else:
            init = sample_ts(s, p["n"], rng, kT=p["kT"], speed=p["speed"], transverse=p["transverse"],
                             hbar=p.get("hbar", 0.0))
        batch.append(dict(surface=s, init=init, gamma=p["gamma"], kT=p["kT"], t_max=p["t_max"],
                          rule=p.get("rule", "first")))
    t0 = time.time()
    res = run_batch(batch, backend=args.backend, seed=args.seed, dt=args.dt)
    out = {"dt": args.dt, "seed": args.seed, "backend": args.backend, "seconds": time.time() - t0, "points": []}
    for p, (o, f) in zip(P, res):
        out["points"].append({**{k: v for k, v in p.items()}, "outcome": ratio(o), "first": ratio(f)})
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(f"{len(P)} points, {sum(p['n'] for p in P)} trajectories, {out['seconds']:.0f} s", flush=True)


def law_microcanonical(s, gamma, E, nq=24):
    """The linear-stage law for fixed energy E on the dividing line (no noise): average over the uniform transverse
    energy E_w of the arcsine law of a fixed-amplitude ellipse, with the speed sqrt(2(E - E_w)) of each slice."""
    Ps, Gs = [], []
    for k in range(nq):
        Ew = E * (k + 0.5) / nq
        L = linear_law(s, float(np.sqrt(2 * (E - Ew))), gamma, 0.0, t_end=40.0 if gamma >= 2 else 12.0)
        om, a, b, m = L["omega_ts"], L["a"], L["b"], L["mean"]
        A = np.sqrt(2 * Ew) * np.sqrt(a * a / om ** 2 + b * b)
        Ps.append(0.5 + np.arcsin(np.clip(m / A, -1, 1)) / np.pi)
        Gs.append(L["G"])
    return {"P_law": float(np.mean(Ps)), "G": float(np.mean(Gs))}


def part_law(args):
    P = points()
    out = {}
    t0 = time.time()
    for p in P:
        s = build(p["surface"])
        t_end = 40.0 if p["gamma"] >= 2 else 12.0
        if p["family"] == "E":
            L = law_microcanonical(s, p["gamma"], p["energy"])
        elif p["family"] == "F":
            L = {**law_dividing_line(s, p["gamma"], p["kT"], t_end=t_end), **tst_split(s, p["kT"])}
        elif p["family"] == "C":
            L = linear_law(s, p["speed"], p["gamma"], 0.0, transverse=p["transverse"], hbar=HBAR, t_end=t_end)
            if p["transverse"] == "qct":
                L = {**L, "P_law": qct_probability(L, HBAR, 0.0)}
        elif p["speed"] is None:
            L = law_thermal_speed(s, p["gamma"], p["kT"], t_end=t_end)
        else:
            L = linear_law(s, p["speed"], p["gamma"], p["kT"], t_end=t_end)
        row = {k: v for k, v in L.items() if k in ("P_law", "z", "G", "mean", "sigma", "sigma_noise", "x_decision",
                                                     "amplification", "P_tst", "y_split")}
        out[p["id"]] = row
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps({"seconds": time.time() - t0, "laws": out}, indent=1))
    print(f"{len(out)} laws, {time.time() - t0:.0f} s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("part", choices=["traj", "law"])
    ap.add_argument("--backend", default="numpy")
    ap.add_argument("--only", default="")
    ap.add_argument("--seed", type=int, default=20261030)
    ap.add_argument("--dt", type=float, default=0.005)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    (part_traj if args.part == "traj" else part_law)(args)


if __name__ == "__main__":
    main()
