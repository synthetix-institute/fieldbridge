"""Figures for the quantum language.

codiscovery_figure(report, path)  (a) the derivation in each realization, letters S A K L; (b) the observable of every
                                  realization that reaches the Bloch rotation, rescaled by its own angle, against the
                                  rotation angle |Omega| t; (c) the dimension of the algebra that the Hamiltonian and the
                                  observable generate
attach_figure(result, path)       the law on the source and on the attached carrier, and the design on that carrier
"""
from __future__ import annotations

from math import cos, pi, sin
from pathlib import Path
from typing import Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

INK, INK2, GRID = "#1f2430", "#6b7280", "#e5e7eb"
BLUE, ORANGE, GREEN, PURPLE, RED, CYAN = "#2563eb", "#ea580c", "#059669", "#7c3aed", "#dc2626", "#0891b2"
LETTER_COLOR = {"S": BLUE, "A": GREEN, "K": PURPLE, "L": CYAN}
COLORS = [BLUE, ORANGE, GREEN, PURPLE, CYAN, RED, INK2]
MARKERS = ["o", "s", "^", "D", "v", "P", "X"]


def _style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7, colors=INK2)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for lab in (ax.xaxis.label, ax.yaxis.label):
        lab.set_fontsize(7.5)
        lab.set_color(INK)


def _title(ax, letter: str, text: str):
    ax.set_title(f"{letter}   {text}", loc="left", fontsize=8.5, fontweight="bold", color=INK)


def _note(step: Dict) -> str:
    L = step["letter"]
    if L == "S":
        return f"dim {step['dimension']} of {step['full_dimension']}"
    if L == "A":
        return f"dim {step['dimension']}"
    if L == "K":
        reps = step["representation"]
        top = max(r["j"] for r in reps)
        return f"j = {int(top) if float(top).is_integer() else f'{int(round(2 * top))}/2'}"
    if L == "L":
        return f"{step['residual']:.0e}"
    return ""


def _derivation_panel(ax, rows):
    from matplotlib.patches import FancyBboxPatch
    heads = {"S": "sector", "A": "algebra", "K": "canonical\nform", "L": "law"}
    xpos = {k: i for i, k in enumerate(heads)}
    for k, x in xpos.items():
        ax.text(x, 0.72, f"{k}\n{heads[k]}", ha="center", va="bottom", fontsize=6.4, color=INK2, linespacing=1.0)
    for i, r in enumerate(rows):
        y = -i
        xs = []
        for st in r["steps"]:
            x = xpos[st["letter"]]
            xs.append(x)
            col = LETTER_COLOR[st["letter"]]
            for fc, ec in ((col, "none"), ("none", col)):
                ax.add_patch(FancyBboxPatch((x - 0.36, y - 0.3), 0.72, 0.6, boxstyle="round,pad=0.02,rounding_size=0.08",
                                            facecolor=fc, alpha=0.14 if fc != "none" else 1.0, edgecolor=ec, lw=0.8))
            ax.text(x, y + 0.08, st["letter"], ha="center", va="center", fontsize=8, fontweight="bold", color=col)
            ax.text(x, y - 0.17, _note(st), ha="center", va="center", fontsize=5.2, color=INK)
        if len(xs) > 1:
            ax.plot(xs, [y] * len(xs), color=GRID, lw=1.2, zorder=0)
        ax.text(-0.6, y + 0.1, r["name"], ha="right", va="center", fontsize=7, color=INK)
        ax.text(-0.6, y - 0.2, r["field"], ha="right", va="center", fontsize=6, color=INK2)
        if r["status"] == "reached":
            ax.text(3.75, y, "reached", ha="left", va="center", fontsize=6.4, color=GREEN)
        else:
            ax.text(max(xs) + 0.55, y, "x", ha="center", va="center", fontsize=8, fontweight="bold", color=RED)
            short = r["obstruction"].split(":")[0].split(";")[0]
            ax.text(3.75, y, short[:70], ha="left", va="center", fontsize=6.2, color=RED)
    ax.set_xlim(-4.2, 8.8)
    ax.set_ylim(-len(rows) + 0.4, 1.45)
    ax.axis("off")


def codiscovery_figure(report: Dict, path: Path):
    rows = report["rows"]
    s = report["summary"]
    reached = [r for r in rows if r["status"] == "reached"]
    ordered = reached + [r for r in rows if r["status"] != "reached"]
    fig = plt.figure(figsize=(12.0, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0])
    a = fig.add_subplot(gs[0, :])
    _derivation_panel(a, ordered)
    _title(a, "a", "Derivation of the Bloch rotation in each realization")
    b = fig.add_subplot(gs[1, 0])
    phi = np.linspace(0, 4 * pi, 300)
    b.plot(phi, np.cos(phi), color=INK, lw=1.0, ls="--", label="cos(|Omega| t)", zorder=3)
    for k, r in enumerate(reached):
        law = r["law"]
        th = r["signature"]["theta_deg"] * pi / 180
        t = np.asarray(law["times"]) * r["signature"]["rate"]
        g = (np.asarray(law["f_exact"]) - cos(th) ** 2) / sin(th) ** 2
        sel = slice(k % 4, None, 4)
        b.plot(t[sel], g[sel], MARKERS[k % len(MARKERS)], ms=3.0, color=COLORS[k % len(COLORS)], lw=0,
               label=f"{r['name']} ({r['field']})")
    b.set_xlabel("rotation angle |Omega| t")
    b.set_ylabel("(f - cos^2 theta) / sin^2 theta,  f = observable / initial value")
    b.legend(fontsize=5.6, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False)
    b.text(0.25, 0.98, f"largest deviation {s['law_residual_max']:.0e}", transform=b.transAxes, fontsize=6.2,
           color=INK2, va="top", ha="center")
    _style(b)
    _title(b, "b", "The observable on every carrier, in the units of the rotation")
    c = fig.add_subplot(gs[1, 1])
    names = [r["name"] for r in ordered][::-1]
    dims = [r["algebra_dimension"] for r in ordered][::-1]
    cols = [GREEN if r["status"] == "reached" else RED for r in ordered][::-1]
    c.barh(range(len(names)), dims, color=cols, alpha=0.75)
    c.axvline(3, color=INK, lw=0.9, ls="--")
    c.text(3.2, len(names) - 0.6, "su(2): dimension 3", fontsize=6.5, color=INK, va="center")
    c.set_yticks(range(len(names)))
    c.set_yticklabels(names, fontsize=6.3)
    c.set_xscale("log", base=2)
    c.set_xticks([1, 2, 3, 4, 8, 16, 24])
    c.set_xticklabels(["1", "2", "3", "4", "8", "16", "24"])
    c.set_xlabel("dimension of the algebra generated by the Hamiltonian and the observable")
    _style(c)
    _title(c, "c", "Algebra generated in each realization")
    fig.suptitle(f"Co-discovery by construction: {s['reached']} realizations from {len(s['fields_reached'])} fields reach "
                 f"the Bloch rotation; {s['obstructed']} are obstructed", fontsize=10, color=INK, x=0.01, ha="left",
                 fontweight="bold")
    fig.savefig(path, dpi=130)
    plt.close(fig)


def attach_figure(result: Dict, path: Path):
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 3.8), constrained_layout=True,
                               gridspec_kw={"width_ratios": [1.4, 1.0]})
    for key, col, mk, ms, face, label in (("source_derivation", BLUE, "o", 6.5, "none", "source"),
                                          ("target_derivation", ORANGE, "s", 2.6, ORANGE, "attached")):
        row = result[key]
        law = row["law"]
        t = np.asarray(law["times"]) * row["signature"]["rate"]
        a.plot(t[::3], np.asarray(law["f_exact"])[::3], mk, ms=ms, color=col, markerfacecolor=face, lw=0,
               label=f"{label}: {row['carrier']}")
    th = result["source_derivation"]["signature"]["theta_deg"] * pi / 180
    phi = np.linspace(0, 4 * pi, 300)
    a.plot(phi, cos(th) ** 2 + sin(th) ** 2 * np.cos(phi), color=INK, lw=1.0, ls="--", label="Rabi law")
    a.set_xlabel("rotation angle |Omega| t")
    a.set_ylabel("observable relative to its initial value")
    a.legend(fontsize=6.2, frameon=False, loc="lower left")
    _style(a)
    _title(a, "a", "The same law on both carriers")
    native = result["native"]
    if "exchange_couplings" in native:
        J = native["exchange_couplings"]
        b.bar(range(len(J)), J, color=ORANGE, alpha=0.8)
        for k, (v, ratio) in enumerate(zip(J, native.get("coupling_ratios", J))):
            b.text(k, v, f"{v:.3f}\n(ratio {ratio:.3f})", ha="center", va="bottom", fontsize=6.5, color=INK)
        b.set_ylim(0, 1.25 * max(J))
        b.set_xticks(range(len(J)))
        b.set_xticklabels([f"{j}-{j + 1}" for j in range(len(J))], fontsize=7)
        b.set_xlabel("bond")
        b.set_ylabel("exchange coupling")
        _style(b)
        _title(b, "b", "Couplings that carry the rotation")
    else:
        b.axis("off")
        b.text(0.02, 0.95, "\n".join(f"{k.replace('_', ' ')}: {v}" for k, v in native.items()), transform=b.transAxes,
               va="top", fontsize=8, color=INK)
        _title(b, "b", "The rotation in the carrier's own operators")
    fig.savefig(path, dpi=130)
    plt.close(fig)
