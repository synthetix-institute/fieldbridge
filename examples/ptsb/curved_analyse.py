"""Figures and summary numbers for the curved valley and the pseudo-minimum (curved_scan.py).

    python curved_analyse.py     # reads results/curved_traj_GH.json, _I.json, _J.json and curved_law_all.json;
                                 # writes figures/curved_*.png and results/curved_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from analyse import CAT, FIG, GRID, HERE, INK, INK2, SEQ  # noqa: E402,F401

T_SLOW = 10.0        # decisions later than this (in units where the passage TS1 -> products takes 1-5) went through
                     # the pseudo-minimum's slow exit


def load():
    pts = {}
    for f in ("curved_traj_GH.json", "curved_traj_I.json", "curved_traj_J.json"):
        d = json.loads((HERE / "results" / f).read_text())
        bins = np.array(d["t_bins"])
        for p in d["points"]:
            p["t_bins"] = bins
            pts[p["id"]] = p
    law = json.loads((HERE / "results/curved_law_all.json").read_text())["laws"]
    for k, p in pts.items():
        row = dict(law[k])
        row["profile_points"] = row.pop("profile", None)      # keep p["profile"], the name of the J profile
        p.update(row)
        o = p["outcome"]
        p["P_out"] = 1 - o["P_plus"]
        p["se"] = o["stderr"]
    return pts


def fam(pts, f, **kw):
    out = [p for p in pts.values() if p["family"] == f]
    for k, v in kw.items():
        out = [p for p in out if (p["surface"].get(k, 0.0) if k in ("kappa", "eps", "a") else p.get(k)) == v]
    return out


MAIN, READ = "law_first_order", "P_law"      # the law shown in the figures; set after the comparison in summary()


def law_of(p, key=None, read=None):
    key = key or MAIN
    if key == "law_reset" and MAIN == "law_first_order":
        key = "law_reset_first_order"
    return p[key].get(read or READ, np.nan) if key in p else np.nan


def slow_share(p):
    """Share of decisions after T_SLOW, and P(outside) among them and among the fast ones."""
    b = p["t_bins"]
    late = b[:-1] >= T_SLOW
    out, inn = np.array(p["time_hist"]["-1"]), np.array(p["time_hist"]["1"])
    n_late, n_fast = out[late].sum() + inn[late].sum(), out[~late].sum() + inn[~late].sum()
    P_late = out[late].sum() / n_late if n_late else np.nan
    P_fast = out[~late].sum() / n_fast if n_fast else np.nan
    return n_late / max(n_late + n_fast, 1), P_late, P_fast, n_late


# ---------------------------------------------------------------------------------------------------- figures
def fig_dynamic_matching(pts):
    kaps = (0.02, 0.05, 0.1, 0.2)
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    ax = axes[0]
    for i, kap in enumerate(kaps):
        rows = sorted(fam(pts, "G", kappa=kap, kT=0.01), key=lambda p: p["gamma"])
        g = np.array([max(p["gamma"], 0.1) for p in rows])
        ax.plot(g, [law_of(p) for p in rows], "-", color=CAT[i], lw=2)
        ax.plot(g, [law_of(p, read="P_wall") for p in rows], "--", color=CAT[i], lw=1)
        ax.errorbar(g, [p["P_out"] for p in rows], yerr=[p["se"] for p in rows], fmt="o", ms=4.5, color=CAT[i],
                    mec="white", mew=0.8, label=f"curvature {kap:g}")
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.3, 1, 3, 10])
    ax.set_xticklabels(["0", "0.3", "1", "3", "10"])
    ax.set_xlabel("friction gamma")
    ax.set_title("kT = 0.01", color=INK, fontsize=9)
    ax.set_ylabel("share of the product outside the bend")
    ax.legend(loc="upper right", fontsize=8)
    for ax, g in zip(axes[1:], (0.0, 1.0)):
        for i, kap in enumerate(kaps):
            rows = sorted(fam(pts, "G", kappa=kap, gamma=g), key=lambda p: p["kT"])
            T = [p["kT"] for p in rows]
            ax.plot(T, [law_of(p) for p in rows], "-", color=CAT[i], lw=2)
            ax.plot(T, [law_of(p, read="P_wall") for p in rows], "--", color=CAT[i], lw=1)
            ax.errorbar(T, [p["P_out"] for p in rows], yerr=[p["se"] for p in rows], fmt="o", ms=4.5,
                        color=CAT[i], mec="white", mew=0.8)
        ax.set_xscale("log")
        ax.set_xlabel("temperature kT")
        ax.set_title(f"friction gamma = {g:g}" + ("  (conservative)" if g == 0 else ""), color=INK, fontsize=9)
    for ax in axes:
        ax.axhline(0.5, color=INK2, lw=0.8, ls=":")
    fig.suptitle("Curved valley with every static quantity symmetric (dotted: IRC, energies, 1/2). Points: "
                 "trajectories; solid: law along the path; dashed: the same law read at the first wall contact",
                 color=INK, fontsize=9.5)
    fig.savefig(FIG / "curved_dynamic_matching.png")
    plt.close(fig)


def fig_against_tilt(pts):
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.4), sharey=True)
    for ax, g in zip(axes, (0.0, 1.0)):
        series = [(0.0, fam(pts, "G", kappa=0.1, gamma=g))] + \
                 [(e, fam(pts, "H", eps=e, gamma=g)) for e in (-0.01, -0.02, -0.04, -0.08)]
        for i, (e, rows) in enumerate(series):
            rows = sorted(rows, key=lambda p: p["kT"])
            T = [p["kT"] for p in rows]
            ax.plot(T, [law_of(p) for p in rows], "-", color=SEQ[i], lw=2)
            ax.errorbar(T, [p["P_out"] for p in rows], yerr=[p["se"] for p in rows], fmt="o", ms=4.5,
                        color=SEQ[i], mec="white", mew=0.8, label=f"tilt inward {abs(e):g}")
        ax.axhline(0.5, color=INK2, lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel("temperature kT")
        ax.set_title(f"friction gamma = {g:g}", color=INK, fontsize=9)
    axes[0].set_ylabel("share of the product outside the bend")
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8)
    fig.suptitle("Curvature 0.1 against a force toward the inside. Points: trajectories; lines: law along the path",
                 color=INK, fontsize=9.5)
    fig.savefig(FIG / "curved_against_tilt.png")
    plt.close(fig)


def fig_pseudo_minimum(pts):
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.2), sharex=True)
    for col, (asym, title) in enumerate((("k", "curved valley, curvature 0.1"), ("e", "straight valley, tilt 0.02"))):
        ax, ax2 = axes[0, col], axes[1, col]
        for i, g in enumerate((0.0, 0.3, 1.0, 3.0)):
            rows = sorted([p for p in fam(pts, "I", gamma=g, kT=0.02) if p["asym"] == asym],
                          key=lambda p: p["surface"]["a"])
            a = [p["surface"]["a"] for p in rows]
            ax.errorbar(a, [p["P_out"] for p in rows], yerr=[p["se"] for p in rows], fmt="o", ms=4.5,
                        color=CAT[i], mec="white", mew=0.8, label=f"gamma = {g:g}")
            ax.plot(a, [law_of(p) for p in rows], "-", color=CAT[i], lw=2)
            ax.plot(a, [law_of(p, "law_reset") for p in rows], "--", color=CAT[i], lw=1.2)
            ax2.plot(a, [slow_share(p)[0] for p in rows], "o-", color=CAT[i], ms=4, lw=1.2)
        ax.axhline(0.5, color=INK2, lw=0.8, ls=":")
        ax.axvline(0.2 - 0.005, color=INK2, lw=0.8, ls="-.")
        ax.set_title(title, color=INK, fontsize=9)
        ax2.set_xlabel("depth of the dip a  (pseudo-minimum for a >= 0.2)")
    axes[0, 0].set_ylabel("share of the product outside")
    axes[1, 0].set_ylabel(f"share decided after t = {T_SLOW:g}")
    axes[0, 1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8)
    fig.suptitle("A dip before the VRI, kT = 0.02. Points: trajectories; solid: law from TS1 (absent where the "
                 "reference stops in the dip); dashed: law from the exit saddle TS1b", color=INK, fontsize=9.5)
    fig.savefig(FIG / "curved_pseudo_minimum.png")
    plt.close(fig)


def fig_held_intermediate(pts):
    profs = (("none", "no intermediate"), ("b4", "exit saddle 0.04 below TS1"), ("b2", "0.02 below"),
             ("b0", "at the level of TS1"))
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.2), sharex=True)
    for col, (asym, title) in enumerate((("k", "curved valley, curvature 0.1"), ("e", "straight valley, tilt 0.02"))):
        ax, ax2 = axes[0, col], axes[1, col]
        for i, (prof, lab) in enumerate(profs):
            rows = sorted([p for p in fam(pts, "J", kT=0.02) if p["asym"] == asym and p["profile"] == prof],
                          key=lambda p: p["gamma"])
            g = np.array([max(p["gamma"], 0.03) for p in rows])
            c = SEQ[1 + i] if i else INK2
            ax.errorbar(g, [p["P_out"] for p in rows], yerr=[p["se"] for p in rows], fmt="o", ms=4.5, color=c,
                        mec="white", mew=0.8, label=lab)
            key = "law" if prof == "none" else "law_reset"
            ax.plot(g, [law_of(p, key) for p in rows], "-" if prof == "none" else "--", color=c, lw=1.6)
            ax2.plot(g, [slow_share(p)[0] for p in rows], "o-", color=c, ms=4, lw=1.2)
        ax.axhline(0.5, color=INK2, lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_xticks([0.03, 0.1, 0.3, 1])
        ax.set_xticklabels(["0", "0.1", "0.3", "1"])
        ax.set_title(title, color=INK, fontsize=9)
        ax2.set_xlabel("friction gamma")
    axes[0, 0].set_ylabel("share of the product outside")
    axes[1, 0].set_ylabel(f"share decided after t = {T_SLOW:g}")
    axes[0, 1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8)
    fig.suptitle("An intermediate 0.1 below TS1, held by an exit saddle near the TS1 level, kT = 0.02. Points: "
                 "trajectories; solid: law from TS1; dashed: law from the exit saddle", color=INK, fontsize=9.5)
    fig.savefig(FIG / "curved_held_intermediate.png")
    plt.close(fig)


# ----------------------------------------------------------------------------------------------------- summary
def summary(pts):
    S = {}
    for f in ("G", "G1", "H", "J"):
        for g in sorted({p["gamma"] for p in fam(pts, f)}):
            rows = fam(pts, f, gamma=g)
            d = {f"{k}:{r}": np.array([p["P_out"] - law_of(p, k, r) for p in rows])
                 for k in ("law", "law_first_order") for r in ("P_law", "P_commit", "P_wall")}
            S[f"{f}_gamma{g:g}"] = {k: {"mean_abs": float(np.nanmean(np.abs(v))), "max_abs": float(np.nanmax(np.abs(v))),
                                        "mean_signed": float(np.nanmean(v))} for k, v in d.items()}
            S[f"{f}_gamma{g:g}"]["n"] = len(rows)
    rows = []
    for p in fam(pts, "I") + fam(pts, "J"):
        sh, P_late, P_fast, n_late = slow_share(p)
        rows.append({"id": p["id"], "family": p["family"], "profile": p.get("profile"), "a": p["surface"].get("a", 0.0), "gamma": p["gamma"], "kT": p["kT"], "asym": p["asym"],
                     "P_out": p["P_out"], "se": p["se"], "law_TS1": law_of(p), "law_TS1b": law_of(p, "law_reset"),
                     "trapped_speeds": p["law"]["trapped"], "exit_barrier": p.get("exit_barrier"),
                     "slow_share": sh, "P_out_slow": P_late, "P_out_fast": P_fast, "n_slow": int(n_late),
                     "reactant": p["outcome"]["reactant"]})
    S["IJ"] = rows
    S["G1"] = [{"id": p["id"], "P_out": p["P_out"], "law": law_of(p), "law_first_order": law_of(p, "law_first_order"),
                "reactant": p["outcome"]["reactant"]} for p in fam(pts, "G1")]
    (HERE / "results/curved_summary.json").write_text(json.dumps(S, indent=1, default=float))
    return S


if __name__ == "__main__":
    pts = load()
    fig_dynamic_matching(pts)
    fig_against_tilt(pts)
    fig_pseudo_minimum(pts)
    fig_held_intermediate(pts)
    S = summary(pts)
    for k, v in S.items():
        if k not in ("IJ", "G1"):
            print(k, v["n"], {kk: {a: round(b, 3) for a, b in vv.items()} for kk, vv in v.items() if kk != "n"})
