"""How much of the product is fixed by the phase-space point at TS1, and by which variables. Cross-validated logistic fits on the starts, against the base rate.

    python analyse_phase_space.py [ens_dir]     # reads ens/decide_E_S.npz from analyse_decision.py; writes results/phase_space.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import log_loss, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from analyse_decision import TIMES  # noqa: E402
from ensemble import Frame  # noqa: E402
from qm import OUT  # noqa: E402


def score(X, y, seed=0):
    """Out-of-fold probabilities over repeated splits: the area under the ROC curve, and the fraction of the base-rate log loss removed."""
    probs = np.zeros((len(y), 4))
    for r in range(4):
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=1, random_state=seed + r)
        for tr, te in cv.split(X, y):
            m = make_pipeline(StandardScaler(), LogisticRegressionCV(Cs=np.logspace(-3, 2, 8), cv=4, max_iter=2000, scoring="neg_log_loss"))
            m.fit(X[tr], y[tr])
            probs[te, r] = m.predict_proba(X[te])[:, 1]
    auc = [roc_auc_score(y, probs[:, r]) for r in range(4)]
    base = log_loss(y, np.full(len(y), y.mean()))
    skill = [1 - log_loss(y, probs[:, r]) / base for r in range(4)]
    return float(np.mean(auc)), float(np.mean(skill))


def main():
    ens = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "ens"
    out = {}
    for e in (31, 32):
        f = Frame(e)
        B = f.C.BL[f.iphi]
        for s in ("cl", "qc"):
            c = np.load(ens / f"decide_{e}_{s}.npz")
            keep = np.isin(c["verdict"], ("Z", "E"))
            y = (c["verdict"][keep] == "Z").astype(int)
            Q, P = c["Q"][keep][:, 1:], c["P"][keep][:, 1:]
            phi = c["phi"][keep]
            rate = P @ B[1:]
            top = np.argsort(-np.abs(B[1:]) / f.w[1:])[:6]                 # the modes that move the torsion most for a thermal amplitude
            sets = {"torsion": phi[:, :1], "torsion+rate": np.c_[phi[:, 0], rate], "torsion+Q0 rate": np.c_[phi[:, 0], rate, c["P"][keep][:, 0]],
                    "soft modes Q,P": np.c_[Q[:, :6], P[:, :6]], "torsion modes Q,P": np.c_[Q[:, top], P[:, top]], "all Q": Q, "all Q,P": np.c_[Q, P]}
            for t in (30.0, 60.0):
                j = int(np.where(TIMES == t)[0][0])
                if np.isfinite(phi[:, j]).all():
                    sets[f"torsion at {t:.0f} fs"] = phi[:, j:j + 1]
            row = dict(n=len(y), z=float(y.mean()), modes=[float(abs(f.wcm[1:][k])) for k in top])
            for k, X in sets.items():
                row[k] = score(X, y)
            out[f"{e}{s}"] = row
            print(f"{e} {s}: n {len(y)} Z {y.mean():.3f}; torsion modes (cm-1) {np.round(row['modes']).astype(int).tolist()}")
            for k in sets:
                print(f"   {k:20s} AUC {row[k][0]:.2f}  skill {row[k][1]:+.2f}", flush=True)
    (OUT / "phase_space.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
