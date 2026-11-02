"""What the thermal ensembles from TS1 show: the product shares, and what at TS1 decides the product.

    python analyse_ensemble.py [ens_dir]        # reads ens/traj/*.done and the starts; writes results/ensemble.json

The product is the torsion well a reactant-to-product trajectory first stays in. Probit fits give the boundary and the width of each start variable:
P(Z | a) = Phi((a* - a) / w), so that a Gaussian spread sigma of a gives the share Phi((a* - <a>) / sqrt(sigma^2 + w^2)).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ensemble import DT_FS, PHI  # noqa: E402
from internal import Curvilinear  # noqa: E402
from qm import BOND, OUT  # noqa: E402
from yamaki_geometries import atoms, coords  # noqa: E402

PUBLISHED = {31: (129, 6), 32: (93, 20)}          # Z and E counts of Wu, Zhang and Hong (2024)


def wilson(k, n, z=1.0):
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (p, c - h, c + h)


def probit(X, y, mult=None):
    """Maximum-likelihood probit of y on X with an intercept. X is standardised inside; returns intercept and slopes in the units of X."""
    X = np.atleast_2d(X.T).T
    mult = np.ones(len(y)) if mult is None else mult
    mu, sd = X.mean(0), X.std(0)
    Z = (X - mu) / sd

    def nll(b):
        s = b[0] + Z @ b[1:]
        return -(mult * (y * norm.logcdf(s) + (1 - y) * norm.logcdf(-s))).sum() + 1e-3 * (b[1:] ** 2).sum()

    r = minimize(nll, np.zeros(X.shape[1] + 1), method="BFGS")
    slopes = r.x[1:] / sd
    return r.x[0] - (slopes * mu).sum(), slopes, r.fun


def phi_dot(C, iphi, q, v):
    """Rate of the torsion at q for the velocity v, in degrees per fs."""
    h = 1e-3
    d = C.wrap(C.xi(q + h * v) - C.xi(q - h * v))[iphi] / (2 * h)
    return np.degrees(d) / 1.03275


def boundary_with_error(phi, y, mult, seed, n=400):
    """The midpoint and width of the probit in the starting torsion, with bootstrap errors."""
    rng = np.random.default_rng(seed)
    a, s, _ = probit(phi, y, mult)
    mid, width = -a / s[0], 1 / abs(s[0])
    boot = []
    p = mult / mult.sum()
    for _ in range(n):
        k = rng.choice(len(y), len(y), p=p)
        if 0 < y[k].sum() < len(y):
            ai, si, _ = probit(phi[k], y[k])
            if si[0] != 0:
                boot.append((-ai / si[0], 1 / abs(si[0])))
    boot = np.array(boot)
    return float(mid), float(width), float(np.median(np.abs(boot[:, 0] - np.median(boot[:, 0]))) * 1.4826), float(np.median(np.abs(boot[:, 1] - np.median(boot[:, 1]))) * 1.4826)


def read_set(ens, e, s):
    """One row per start that has a forward verdict."""
    ini = np.load(ens / f"init_{e}_{s}.npz")
    if s == "cl" and (ens / f"order_{e}_cl.npz").exists():
        order = np.load(ens / f"order_{e}_cl.npz")["order"]
    else:
        order = np.arange(len(ini["q"]))
    rows = {}
    fails = 0
    # a repeated draw counts only up to the first start that has not finished: later repeats would favour the heavy starts that happen to have run
    for draw, i in enumerate(order):
        name = f"{e}{s}_{int(i):04d}"
        if i in rows:
            rows[i]["mult"] += 1
            continue
        fd, bd = ens / "traj" / f"{name}_f.done", ens / "traj" / f"{name}_b.done"
        if (ens / "traj" / f"{name}_f.failed").exists():
            fails += 1
            continue
        if not fd.exists():
            break
        rows[i] = dict(i=int(i), draw=draw, mult=1, forward=fd.read_text(), backward=bd.read_text() if bd.exists() else None)
    return ini, rows, fails


def main():
    ens = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "ens"
    out = {}
    for e in (31, 32):
        ts = np.load(OUT / f"ts_{e}.npz")
        irc = np.load(OUT / f"irc_{e}.npz")
        C = Curvilinear([a[0] for a in atoms(e, "TS1")], [coords(e, n) for n in ("TS1", "PM", "P1", "P2")], BOND[e], irc["q"][0], ts["m"], ts["L"])
        iphi = C.nb + C.na + C.g[2].index(PHI[e])
        phi_ts = float(np.degrees(C.xi(irc["q"][0])[iphi]))
        for s in ("cl", "qc"):
            ini, rows, fails = read_set(ens, e, s)
            if not rows:
                continue
            rows = {i: r for i, r in rows.items() if r["backward"] is not None or r["forward"] == "reactant"}   # complete pairs only
            idx = np.array(sorted(rows))
            verdict = [rows[i]["forward"] for i in idx]
            back = [rows[i]["backward"] for i in idx]
            mult = np.array([rows[i]["mult"] for i in idx], float)
            phi = np.degrees(ini["phi"][idx])
            phi = phi_ts + (phi - phi_ts + 180) % 360 - 180          # continuous through +-180 on the Z side, relative to the torsion at TS1
            pdot = np.array([phi_dot(C, iphi, ini["q"][i], ini["v"][i]) for i in idx])
            count = {k: float(mult[[v == k for v in verdict]].sum()) for k in ("Z", "E", "reactant", "undecided")}
            # reactant-to-product: the backward segment must return to the reactant side
            reactive = np.array([v in ("Z", "E") and b != "product" for v, b in zip(verdict, back)])
            z = np.array([v == "Z" for v in verdict])
            nz, ne = mult[reactive & z].sum(), mult[reactive & ~z].sum()
            share = wilson(nz, nz + ne)
            row = dict(entry=e, sampling=s, starts=len(idx), draws=float(mult.sum()), failed=fails, counts=count,
                       backward_done=int(sum(b is not None for b in back)), both_product=float(mult[[b == "product" and v in ("Z", "E") for v, b in zip(verdict, back)]].sum()),
                       z=float(nz), e=float(ne), share=list(map(float, share)), phi_ts=phi_ts,
                       phi_mean=float(np.average(phi[reactive], weights=mult[reactive])), phi_sd=float(np.sqrt(np.cov(phi[reactive], aweights=mult[reactive]))))
            if nz >= 3 and ne >= 3:
                y = z[reactive].astype(float)
                m = mult[reactive]
                mid, width, emid, ewidth = boundary_with_error(phi[reactive], y, m, e)
                row.update(boundary=mid, boundary_err=emid, width=width, width_err=ewidth)
                row["share_from_gaussian"] = float(norm.cdf((mid - row["phi_mean"]) / np.hypot(row["phi_sd"], width)))
                row["share_below_boundary"] = float(m[phi[reactive] < mid].sum() / m.sum())
                a2, s2, _ = probit(np.c_[phi[reactive], pdot[reactive]], y, m)
                row["torsion_and_rate"] = dict(boundary=float(-a2 / s2[0]), lookahead_fs=float(s2[1] / s2[0]), width=float(1 / abs(s2[0])), rate_sd=float(pdot[reactive].std()))
                ok = phi[reactive] < mid
                row["accuracy_boundary"] = float(m[ok == (y > 0.5)].sum() / m.sum())
            # bonding time and first choice time from the stored frames, and the energy drift
            drift, times = [], []
            for i in idx:
                d = np.load(ens / "traj" / f"{e}{s}_{int(i):04d}_f.npz")
                tot = d["E"] + 0.5 * (d["v"] ** 2).sum(1)
                drift.append(float(np.abs(tot - tot[0]).max() * 627.509))
                times.append(float((len(d["E"]) - 1) * DT_FS))
            row["energy_drift_kcal_median"] = float(np.median(drift))
            row["energy_drift_kcal_max"] = float(np.max(drift))
            row["forward_fs_mean"] = float(np.mean(times))
            row["published_share"] = PUBLISHED[e][0] / sum(PUBLISHED[e])
            row["table"] = dict(i=idx.tolist(), mult=mult.tolist(), verdict=verdict, backward=back, phi=phi.round(2).tolist(), phidot=pdot.round(4).tolist(), reactive=reactive.tolist())
            out[f"{e}{s}"] = row
            print(f"{e} {s}: {len(idx)} starts, {fails} failed; Z {count['Z']:.0f} E {count['E']:.0f} back {count['reactant']:.0f} undecided {count['undecided']:.0f}; backward done {row['backward_done']}, product at both ends {row['both_product']:.0f}; "
                  f"share Z {share[0]:.3f} ({share[1]:.3f}-{share[2]:.3f}), published {row['published_share']:.3f}; phi {row['phi_mean']:.1f} +- {row['phi_sd']:.1f}; drift median {row['energy_drift_kcal_median']:.2f} max {row['energy_drift_kcal_max']:.2f}; forward {row['forward_fs_mean']:.0f} fs")
            if "boundary" in row:
                print(f"   boundary {row['boundary']:.1f} +- {row['boundary_err']:.1f}, width {row['width']:.1f} +- {row['width_err']:.1f}; Gaussian share {row['share_from_gaussian']:.3f}; below boundary {row['share_below_boundary']:.3f}; accuracy {row['accuracy_boundary']:.3f}; with rate: {row['torsion_and_rate']}")
    (OUT / "ensemble.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
