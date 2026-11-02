"""Whether a product reached first is left again: the torsion of the continued trajectories up to 400 fs.

    python analyse_late.py [ens_dir]        # reads ens/late/*.npz; writes results/late.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ensemble import DT_FS  # noqa: E402
from qm import OUT  # noqa: E402

HOLD_FS = 40.0


def wells(phi, R):
    """The sequence of wells, each counted when the torsion stays within 73 degrees of planar on that side for HOLD_FS after coming within 45 of it, as in the forward rule."""
    c = np.cos(phi)
    seq, side, since = [], None, 0
    for k in range(len(phi)):
        s = side
        if s == "Z" and c[k] > -0.3 or s == "E" and c[k] < 0.3 or R[k] > 1.8:
            s = None
        if s is None and R[k] < 1.8:
            s = "Z" if c[k] < -0.7 else "E" if c[k] > 0.7 else None
        if s != side:
            side, since = s, k
        if side is not None and (k - since) * DT_FS >= HOLD_FS and (not seq or seq[-1] != side):
            seq.append(side)
    return seq


def main():
    ens = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "ens"
    rows = []
    for p in sorted((ens / "late").glob("late_*.done")):
        d = np.load(p.with_suffix(".npz"))
        name = p.stem.split("_")
        phi, R, E = d["phi"], d["R"], d["E"]
        seq = wells(phi, R)
        first = str(d["first"])
        rows.append(dict(set=name[1], i=int(name[2]), first=first, start_fs=float(d["start_fs"]), wells=seq, end=seq[-1] if seq else first, left=any(s != first for s in seq), R_max=float(R.max()), drift_kcal=float(np.abs(E - E[0]).max() * 627.509), rms_fs=float((len(phi)) * DT_FS)))
    out = {}
    for s in sorted({r["set"] for r in rows}):
        for side in ("Z", "E"):
            sub = [r for r in rows if r["set"] == s and r["first"] == side]
            left = sum(r["left"] for r in sub)
            end = sum(r["end"] != side for r in sub)
            fs = sum(r["rms_fs"] for r in sub)
            out[f"{s}{side}"] = dict(n=len(sub), left=left, ended_other=end, fs=fs)
            print(f"{s} first {side}: {len(sub)} trajectories, {fs:.0f} fs in all, left the first well in {left}, ended in the other in {end}; largest bond {max(r['R_max'] for r in sub):.2f}; drift max {max(r['drift_kcal'] for r in sub):.2f}")
    out["rows"] = rows
    (OUT / "late.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
