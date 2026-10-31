"""When the product is decided: how well the torsion at each time after TS1 separates Z from E, and whether free motion of the TS1 modes accounts for it.

    python analyse_decision.py [ens_dir]        # reads ens/traj/*_f.npz and the starts; caches ens/decide_E_S.npz; writes results/decision.json

The torsion at time t is compared with the torsion the TS1 modes would give by harmonic motion alone,
    phi_h(t) = phi(0) + sum_k B_k [ Q_k (cos w_k t - 1) + P_k sin(w_k t) / w_k ],
with B_k the derivative of the torsion along mode k at TS1, and with the reaction mode left out of the sum.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from analyse_ensemble import probit  # noqa: E402
from ensemble import DT_FS, TIME_UNIT_FS, Frame  # noqa: E402
from qm import OUT, project  # noqa: E402

TIMES = np.arange(0, 125, 5.0)


def table(f, ens, e, s):
    """Torsion and bond at the times TIMES for every finished forward trajectory, and the modes and their rates at the start. New trajectories are added to the cache."""
    cache = ens / f"decide_{e}_{s}.npz"
    old = {}
    if cache.exists():
        c = np.load(cache, allow_pickle=False)
        old = {int(i): k for k, i in enumerate(c["i"])}
    phi_ts = np.degrees(f.phi(f.q0))
    rows = []
    for p in sorted((ens / "traj").glob(f"{e}{s}_*_f.done")):
        i = int(p.name.split("_")[1])
        if i in old:
            k = old[i]
            rows.append((i, c["verdict"][k], c["phi"][k], c["R"][k], c["Q"][k], c["P"][k]))
            continue
        d = np.load(ens / "traj" / f"{e}{s}_{i:04d}_f.npz")
        n = len(d["q"])
        phi, R = np.full(len(TIMES), np.nan), np.full(len(TIMES), np.nan)
        steps = [int(round(t / DT_FS)) for t in TIMES]
        last = None
        for j, k in enumerate(steps):
            if k >= n:
                break
            a = np.degrees(f.phi(d["q"][k]))
            ref = phi_ts if last is None else last
            a = ref + (a - ref + 180) % 360 - 180             # continuous in time, relative to TS1 at the start
            phi[j], R[j], last = a, f.bond(d["q"][k]), a
        q0, v0 = d["q"][0], d["v"][0]
        A = f.C.pinv @ f.bmatrix(q0) @ project(q0, f.m)
        rows.append((i, p.read_text(), phi, R, f.C(q0), A @ v0))
    if rows:
        np.savez(cache, i=np.array([r[0] for r in rows]), verdict=np.array([r[1] for r in rows]), phi=np.array([r[2] for r in rows]), R=np.array([r[3] for r in rows]),
                 Q=np.array([r[4] for r in rows]), P=np.array([r[5] for r in rows]))
    return rows


def harmonic(f, iphi, phi0, Q, P, t_fs):
    """phi_h(t), in degrees, without the reaction mode."""
    t = t_fs / TIME_UNIT_FS
    B = f.C.BL[iphi]
    w = f.w
    out = np.radians(phi0).copy()
    for k in range(1, len(w)):
        out = out + B[k] * (Q[:, k] * (np.cos(w[k] * t) - 1) + P[:, k] * np.sin(w[k] * t) / w[k])
    return np.degrees(out)


def main():
    ens = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "ens"
    out = {}
    for e in (31, 32):
        f = Frame(e)
        for s in ("cl", "qc"):
            rows = table(f, ens, e, s)
            if not rows:
                continue
            v = np.array([r[1] for r in rows])
            phi, R, Q, P = (np.array([r[k] for r in rows]) for k in (2, 3, 4, 5))
            keep = np.isin(v, ("Z", "E"))
            y = (v[keep] == "Z").astype(float)
            curve = []
            for j, t in enumerate(TIMES):
                ok = np.isfinite(phi[keep, j])
                if ok.sum() < 8 or not 0 < y[ok].sum() < ok.sum():
                    continue
                a, b, _ = probit(phi[keep, j][ok], y[ok])
                h = harmonic(f, f.iphi, phi[keep, 0][ok], Q[keep][ok], P[keep][ok], t)
                best = max(((phi[keep, j][ok] < c) == (y[ok] > 0)).mean() for c in phi[keep, j][ok])
                besth = max(((h < c) == (y[ok] > 0)).mean() for c in h)
                r = np.corrcoef(h, phi[keep, j][ok])[0, 1] if t > 0 else 1.0
                curve.append(dict(t=float(t), n=int(ok.sum()), R=float(np.nanmean(R[keep, j][ok])), mean=float(phi[keep, j][ok].mean()), sd=float(phi[keep, j][ok].std()),
                                  mid=float(-a / b[0]), width=float(1 / abs(b[0])), accuracy=float(best), accuracy_harmonic=float(besth), r_harmonic=float(r), sd_harmonic=float(h.std())))
            out[f"{e}{s}"] = dict(n=int(keep.sum()), z=float(y.sum()), curve=curve)
            print(f"{e} {s}: {keep.sum()} products, Z {y.mean():.3f}")
            for c in curve:
                if c["t"] % 10 == 0:
                    print(f"  {c['t']:5.0f} fs  n {c['n']:3d}  R {c['R']:.2f}  torsion {c['mean']:7.1f} +- {c['sd']:5.1f}  boundary {c['mid']:7.1f} width {c['width']:6.1f}  accuracy {c['accuracy']:.2f}; free modes: accuracy {c['accuracy_harmonic']:.2f}, corr {c['r_harmonic']:.2f}, sd {c['sd_harmonic']:.1f}")
    (OUT / "decision.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
