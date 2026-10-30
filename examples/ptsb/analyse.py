"""Figures and summary numbers from results/traj.json, results/law.json and results/classes.json.

    python analyse.py            # writes figures/*.png and results/summary.json
"""
from __future__ import annotations

import json
from pathlib import Path
from statistics import NormalDist

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
SEQ = ["#9ec5f4", "#6da7ec", "#3987e5", "#1c5cab", "#0d366b"]          # one hue, light to dark: temperature
inv = NormalDist().inv_cdf

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})


def load():
    traj = {p["id"]: p for p in json.loads((HERE / "results/traj.json").read_text())["points"]}
    law = json.loads((HERE / "results/law.json").read_text())["laws"]
    for k, p in traj.items():
        p["law"] = law[k]
    return traj


def probit(P):
    return inv(min(max(P, 1e-4), 1 - 1e-4))


def fam(traj, f):
    return [p for p in traj.values() if p["family"] == f]


def fig_asymmetry_temperature(traj):
    A = fam(traj, "A")
    gammas = (0.0, 1.0, 10.0)
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.3), sharey=True)
    for ax, g in zip(axes, gammas):
        for i, kT in enumerate((0.005, 0.01, 0.02, 0.05, 0.1)):
            rows = sorted((p for p in A if p["gamma"] == g and p["kT"] == kT), key=lambda p: p["surface"][1])
            e = np.array([p["surface"][1] for p in rows])
            P = 1 - np.array([p["first"]["P_plus"] for p in rows])          # share of the favoured product
            se = np.array([p["first"]["stderr"] for p in rows])
            L = 1 - np.array([p["law"]["P_law"] for p in rows])
            ax.plot(e, L, "-", color=SEQ[i], lw=2)
            ax.errorbar(e, P, yerr=se, fmt="o", ms=4.5, color=SEQ[i], mec="white", mew=0.8, label=f"kT = {kT:g}")
        ax.set_xscale("log")
        ax.set_title(f"friction gamma = {g:g}" + ("  (conservative)" if g == 0 else ""), color=INK, fontsize=9)
        ax.set_xlabel("tilt across the path, eps")
    axes[0].set_ylabel("share of the favoured product")
    axes[0].legend(loc="upper left", fontsize=8)
    fig.suptitle("Minimal valley-ridge surface: points, trajectories (first entry); lines, write law along the path",
                 color=INK, fontsize=10)
    fig.savefig(FIG / "ratio_vs_asymmetry_temperature.png")
    plt.close(fig)


def fig_law_vs_traj(traj):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.6))
    out = {}
    for fam_id, marker, name in (("A", "o", "minimal surface"), ("D", "s", "Collins et al. surface")):
        rows = [p for p in fam(traj, fam_id) if np.isfinite(p["law"]["P_law"])]
        gs = sorted({p["gamma"] for p in rows})
        ratios = {}
        for i, g in enumerate(gs):
            r = [p for p in rows if p["gamma"] == g]
            zl = np.array([probit(p["law"]["P_law"]) for p in r])
            zt = np.array([probit(p["first"]["P_plus"]) for p in r])
            col = CAT[min(i, len(CAT) - 1)] if fam_id == "A" else INK2
            if fam_id == "A":
                a1.plot(zl, zt, marker, ms=3.5, color=col, mec="white", mew=0.5, label=f"gamma = {g:g}")
            else:
                a1.plot(zl, zt, marker, ms=4, mfc="none", color=col, mew=0.8,
                        label="Collins et al., all gamma" if i == 0 else None)
            ok = np.abs(zl) > 0.1
            ratios[g] = (float(np.median(zt[ok] / zl[ok])), float(np.percentile(zt[ok] / zl[ok], 25)),
                         float(np.percentile(zt[ok] / zl[ok], 75)), int(ok.sum()))
        out[fam_id] = ratios
        g = np.array(list(ratios))
        med = np.array([ratios[k][0] for k in ratios])
        lo = np.array([ratios[k][1] for k in ratios])
        hi = np.array([ratios[k][2] for k in ratios])
        gp = np.where(g == 0, 0.03, g)
        a2.errorbar(gp, med, yerr=[med - lo, hi - med], fmt=marker + "-", ms=5, lw=1.5,
                    color=CAT[0] if fam_id == "A" else CAT[1], mfc="white" if fam_id == "D" else None, label=name)
    lim = 1.7
    a1.plot([-lim, lim], [-lim, lim], "-", color=INK2, lw=0.8)
    a1.set_xlim(-lim, lim)
    a1.set_ylim(-lim, lim)
    a1.set_xlabel("probit of the law, z_law")
    a1.set_ylabel("probit of the trajectories, z_traj")
    a1.legend(fontsize=7, loc="upper left")
    a1.set_title("law against trajectories (first entry)", color=INK, fontsize=9)
    a2.axhline(1.0, color=INK2, lw=0.8)
    a2.set_xscale("log")
    a2.set_xlabel("friction gamma (0 plotted at 0.03)")
    a2.set_ylabel("z_traj / z_law (median, quartiles)")
    a2.legend(fontsize=8)
    a2.set_title("excess selectivity of the conservative passage", color=INK, fontsize=9)
    fig.savefig(FIG / "law_vs_trajectories.png")
    plt.close(fig)
    return out


def fig_zero_point(traj):
    C = fam(traj, "C")
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.3), sharey=True)
    for ax, g in zip(axes, (0.0, 3.0)):
        for i, (tr, name) in enumerate((("qct", "normal-mode (fixed energy, uniform phase)"),
                                        ("wigner", "Gaussian (Wigner), same mean energy"))):
            rows = sorted((p for p in C if p["gamma"] == g and p["transverse"] == tr), key=lambda p: p["surface"][1])
            e = np.array([p["surface"][1] for p in rows])
            ax.plot(e, 1 - np.array([p["law"]["P_law"] for p in rows]), "-", color=CAT[i], lw=2)
            ax.errorbar(e, 1 - np.array([p["first"]["P_plus"] for p in rows]),
                        yerr=[p["first"]["stderr"] for p in rows], fmt="o", ms=4.5, color=CAT[i], mec="white",
                        mew=0.8, label=name)
        ax.set_xscale("log")
        ax.set_xlabel("tilt across the path, eps")
        ax.set_title(f"zero-point motion only, gamma = {g:g}", color=INK, fontsize=9)
    axes[0].set_ylabel("share of the favoured product")
    axes[0].legend(fontsize=7, loc="upper left")
    fig.suptitle("points: trajectories; lines: arcsine law (normal-mode) and Gaussian law (Wigner) along the path",
                 color=INK, fontsize=9)
    fig.savefig(FIG / "zero_point_sampling.png")
    plt.close(fig)


def fig_collins(traj):
    E = fam(traj, "E")
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
    for ax, En in zip(axes, (0.01, 0.1)):
        for i, xs in enumerate((2.0, 2.05, 2.1)):
            rows = sorted((p for p in E if p["energy"] == En and p["surface"][1] == xs), key=lambda p: p["gamma"])
            g = np.array([p["gamma"] for p in rows])
            ax.plot(g, [p["outcome"]["P_plus"] for p in rows], "o-", ms=4, lw=1.5, color=CAT[i],
                    label=f"x* = {xs:g}, final well")
            ax.plot(g, [p["first"]["P_plus"] for p in rows], "^", ms=4, mfc="white", color=CAT[i],
                    label=f"x* = {xs:g}, first entry")
            ax.plot(g, [p["law"]["P_law"] for p in rows], "--", lw=1.2, color=CAT[i])
        ax.set_xlabel("damping gamma (no noise)")
        ax.set_title(f"Collins et al. protocol, E = {En:g} above TS1", color=INK, fontsize=9)
    axes[0].set_ylabel("share of the upper product")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=7, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.12))
    fig.suptitle("filled: well occupied at the end (their definition); open: first entry; dashed: write law",
                 color=INK, fontsize=9)
    fig.savefig(FIG / "collins_protocol.png")
    plt.close(fig)


def fig_cusp(traj):
    F = fam(traj, "F")
    fig, ax = plt.subplots(figsize=(6, 3.6))
    for i, kT in enumerate((0.01, 0.03, 0.1)):
        col = SEQ[[1, 3, 4][i]]
        for g, mk in ((0.0, "o"), (3.0, "s")):
            rows = sorted((p for p in F if p["kT"] == kT and p["gamma"] == g), key=lambda p: p["surface"][3])
            d = np.array([p["surface"][3] for p in rows])
            ax.plot(d, [1 - p["first"]["P_plus"] for p in rows], mk, ms=4.5, color=col,
                    mfc=col if g == 0 else "white", mew=0.9,
                    label=f"kT = {kT:g}, gamma = {g:g}")
        rows = sorted((p for p in F if p["kT"] == kT and p["gamma"] == 0.0), key=lambda p: p["surface"][3])
        d = np.array([p["surface"][3] for p in rows])
        law = np.array([1 - p["law"]["P_law"] for p in rows])
        tst = np.array([1 - p["law"]["P_tst"] for p in rows])
        ax.plot(d, law, "-", color=col, lw=1.8)
        ax.plot(d, tst, ":", color=col, lw=1.8)
    ax.axvline(-0.25, color=INK2, lw=0.8, ls="--")
    ax.text(-0.247, 0.98, "VRI on TS line (eps = 0)", rotation=90, ha="left", va="top", color=INK2, fontsize=7,
            transform=ax.get_xaxis_transform())
    ax.axvline(-0.2865, color=INK2, lw=0.8)
    ax.text(-0.289, 0.98, "fold on TS line", rotation=90, ha="right", va="top", color=INK2,
            fontsize=7, transform=ax.get_xaxis_transform())
    ax.text(-0.44, 0.52, "two saddles\n(small barrier)", color=INK2, fontsize=8, va="top")
    ax.text(0.0, 0.52, "one saddle, VRI after it\n(downhill)", color=INK2, fontsize=8, va="top")
    ax.set_xlabel("position of the VRI relative to the TS line, d (VRI at x = 1/4 + d)")
    ax.set_ylabel("share of the favoured product")
    ax.set_title("eps = 0.02: solid, write law; dotted, transition-state split", color=INK, fontsize=9)
    ax.legend(fontsize=6.5, loc="upper right", ncol=2)
    ax.set_ylim(0.45, 1.0)
    fig.savefig(FIG / "cusp_crossover.png")
    plt.close(fig)


def fig_cubic(traj):
    B = fam(traj, "B")
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for i, c in enumerate((0.1, 0.2, 0.4)):
        for g, mk in ((0.0, "o"), (3.0, "s")):
            rows = sorted((p for p in B if p["surface"][2] == c and p["gamma"] == g), key=lambda p: p["kT"])
            ax.errorbar([p["kT"] for p in rows], [p["first"]["P_plus"] for p in rows],
                        yerr=[p["first"]["stderr"] for p in rows], fmt=mk + "-", ms=4, lw=1, color=CAT[i],
                        mfc=CAT[i] if g == 0 else "white", label=f"c = {c:g}, gamma = {g:g}")
    ax.axhline(0.5, color=INK2, lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("kT")
    ax.set_ylabel("share of the shallower (upper) product")
    ax.set_title("asymmetry of the valleys, no force on the path", color=INK, fontsize=9)
    ax.legend(fontsize=6.5, ncol=2, loc="lower left")
    ax.set_ylim(0.4, 0.56)
    fig.savefig(FIG / "valley_shape_asymmetry.png")
    plt.close(fig)


def summary(traj, ratios):
    out = {"z_ratio_by_gamma": ratios}
    for f in ("A", "A1", "B", "C", "D", "E", "F"):
        rows = [p for p in fam(traj, f) if np.isfinite(p["law"].get("P_law", np.nan))]
        dP = np.array([p["first"]["P_plus"] - p["law"]["P_law"] for p in rows])
        out[f] = {"points": len(rows), "mean_abs_dP": float(np.mean(np.abs(dP))) if rows else None,
                  "max_abs_dP": float(np.max(np.abs(dP))) if rows else None,
                  "reactant_share_mean": float(np.mean([p["first"]["reactant"] / p["first"]["n"] for p in rows]))
                  if rows else None}
        for g in sorted({p["gamma"] for p in rows}):
            r = [p for p in rows if p["gamma"] == g]
            d = np.array([p["first"]["P_plus"] - p["law"]["P_law"] for p in r])
            out[f][f"gamma_{g:g}"] = {"points": len(r), "mean_abs_dP": float(np.mean(np.abs(d))),
                                      "max_abs_dP": float(np.max(np.abs(d)))}
    F = fam(traj, "F")
    tst = [p for p in F if np.isfinite(p["law"]["P_tst"])]
    out["F_tst"] = {"points": len(tst), "mean_abs_dP": float(np.mean([abs(p["first"]["P_plus"] - p["law"]["P_tst"])
                                                                      for p in tst]))}
    (HERE / "results/summary.json").write_text(json.dumps(out, indent=1))
    return out


def main():
    FIG.mkdir(exist_ok=True)
    traj = load()
    fig_asymmetry_temperature(traj)
    ratios = fig_law_vs_traj(traj)
    fig_zero_point(traj)
    fig_collins(traj)
    fig_cusp(traj)
    fig_cubic(traj)
    print(json.dumps(summary(traj, ratios), indent=1))


if __name__ == "__main__":
    main()
