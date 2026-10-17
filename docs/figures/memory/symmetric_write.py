"""Manuscript figure: one write mechanism, the symmetric write, in eight physical systems from six fields.

  a  each system with the parameter that is swept and the two states it writes (blue: the state the bias favours);
     green, the two systems that reach the mechanism only after a second parameter is tuned to a cusp
  b  the landscape of each system at its write point, V = -int g dx in the canonical coordinate x = |a3|^(1/2) s,
     against the quartic x^4/4 of the pitchfork; insets: the landscape before, at and after the write point
  c  a single-mode laser swept through threshold with a weak seed field (the protocol of the law test at a predicted
     accuracy of 0.8, Gamma = 0.05, r = 0.02): the field stays near zero past threshold and then grows with phase 0
     or pi; black, the static branches; grey band, the window of the linear stage in which the bias acts, +- r^(1/2)
     in the growth rate of the critical mode
  d  the accuracy of the write against the bias theta = h_s / (D_s^(1/2) r^(1/4)) for all eight systems (six biases,
     2000 trajectories each, the protocol of construct.swept_write_check with the signs of the derivation), against
     P = Phi(pi^(1/4) theta); inset: the constant of each system from docs/site/law_constants.json (25,600
     trajectories)

usage (from the repository root, with python3.11):
  python docs/figures/memory/symmetric_write.py compute   # about 15 minutes; writes symmetric_write_data.json
  python docs/figures/memory/symmetric_write.py draw      # writes symmetric_write.pdf, .png and .svg
"""
import json
import sys
import time
from math import erf, pi, sqrt
from pathlib import Path
from statistics import NormalDist

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches as mp
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
SEED = 20260923
DATA = HERE / "symmetric_write_data.json"
STEM = HERE / "symmetric_write"

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, RULE, FILL = "#0b0b0b", "#52514e", "#8a8984", "#d9d8d3", "#f0efec"
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 7, "mathtext.fontset": "stixsans",
    "axes.linewidth": 0.6, "axes.edgecolor": INK2, "axes.labelcolor": INK, "axes.labelsize": 7,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "text.color": INK, "pdf.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False})

LAW = pi ** 0.25
MODELS = ["pitchfork", "laser", "convection_roll", "stoner_wohlfarth", "toggle", "repressor_ring4", "schlogl",
          "toggle_unequal"]
META = {
    "pitchfork": dict(name="Landau theory", field="statistical physics", sweep=r"temperature, $\varepsilon$", marker="o",
                      states=("s > 0", "s < 0")),
    "laser": dict(name="Laser", field="laser physics", sweep="pump $P$", marker="s", states=("phase 0", "phase π")),
    "convection_roll": dict(name="Convection roll", field="fluid dynamics", sweep="heating $r$", marker="D",
                            states=("cw roll", "ccw roll")),
    "stoner_wohlfarth": dict(name="Nanomagnet", field="magnetism", sweep="field $h$", marker="^",
                             states=("up", "down")),
    "toggle": dict(name="Toggle switch", field="synthetic biology", sweep=r"expression $\alpha$", marker="v",
                   states=("A on", "B on")),
    "repressor_ring4": dict(name="Repressor ring", field="synthetic biology", sweep=r"expression $\alpha$", marker="p",
                            states=("1010", "0101")),
    "schlogl": dict(name="Schlögl reactor", field="chemical kinetics", sweep="feed $b$, at a cusp", marker="h",
                    states=("high X", "low X"), cusp=True),
    "toggle_unequal": dict(name="Unequal toggle", field="synthetic biology", sweep=r"$\alpha$, at a cusp", marker="X",
                           states=("A on", "B on"), cusp=True),
}


def color_of(key):
    return AQUA if META[key].get("cusp") else INK2


def phi(x):
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))


# ------------------------------------------------------------------------------------------------ pictograms
def arrow(ax, a, b, color=INK2, lw=0.8, style="-|>", ms=5, **kw):
    ax.add_patch(mp.FancyArrowPatch(a, b, arrowstyle=style, mutation_scale=ms, color=color, lw=lw,
                                    shrinkA=0, shrinkB=0, **kw))


def repress(ax, a, b, rad, color=INK2, lw=0.8):
    ax.add_patch(mp.FancyArrowPatch(a, b, arrowstyle=mp.ArrowStyle("-[", widthB=3.0, lengthB=0.0), lw=lw,
                                    color=color, connectionstyle=f"arc3,rad={rad}", shrinkA=1, shrinkB=1))


def ellipse_arc(ax, c, rx, ry, t0, t1, color, lw=0.9, head=True, ms=5):
    t = np.linspace(np.radians(t0), np.radians(t1), 60)
    x, y = c[0] + rx * np.cos(t), c[1] + ry * np.sin(t)
    ax.plot(x[:-3], y[:-3], color=color, lw=lw, solid_capstyle="round")
    if head:
        arrow(ax, (x[-4], y[-4]), (x[-1], y[-1]), color=color, lw=lw, ms=ms)


def draw_landau(ax):
    x = np.linspace(-1.25, 1.25, 120)
    V = (x ** 2 - 0.8) ** 2
    ax.plot(0.5 + 0.33 * x, 0.5 + 0.36 * V / V.max(), color=INK2, lw=1.0)
    ax.add_patch(mp.Circle((0.5, 0.5 + 0.36 * 0.64 / V.max() + 0.045), 0.04, color=MUTED, zorder=3))


def draw_laser(ax):
    for x0 in (0.08, 0.885):
        ax.add_patch(mp.Rectangle((x0, 0.56), 0.035, 0.34, color=INK2, lw=0))
    ax.add_patch(mp.FancyBboxPatch((0.3, 0.66), 0.4, 0.14, boxstyle="round,pad=0.0,rounding_size=0.03",
                                   fc=FILL, ec=INK2, lw=0.6))
    x = np.linspace(0.115, 0.885, 200)
    ax.plot(x, 0.73 + 0.07 * np.sin((x - 0.115) / 0.77 * 5 * pi), color=INK2, lw=0.7)
    for xa in (0.38, 0.5, 0.62):
        arrow(ax, (xa, 0.5), (xa, 0.645), color=MUTED, lw=0.7, ms=4)
    ax.text(0.27, 0.53, "$P$", ha="center", va="center", fontsize=6.5, color=MUTED)


def draw_roll(ax):
    ax.add_patch(mp.Rectangle((0.08, 0.53), 0.84, 0.39, fc="none", ec=INK2, lw=0.6))
    ax.plot([0.08, 0.92], [0.53, 0.53], color=INK2, lw=1.8)
    for rx, ry in ((0.3, 0.14), (0.17, 0.075)):
        ax.add_patch(mp.Ellipse((0.5, 0.725), 2 * rx, 2 * ry, fc="none", ec=MUTED, lw=0.8))
    for xh in (0.25, 0.5, 0.75):
        yy = np.linspace(0.45, 0.515, 20)
        ax.plot(xh + 0.012 * np.sin((yy - 0.45) / 0.065 * 2 * pi), yy, color=MUTED, lw=0.6)


def draw_magnet(ax):
    ax.add_patch(mp.Ellipse((0.58, 0.72), 0.24, 0.42, fc=FILL, ec=INK2, lw=0.7))
    ax.plot([0.58, 0.58], [0.48, 0.96], color=MUTED, lw=0.6, ls=(0, (2, 2)))
    arrow(ax, (0.58, 0.72), (0.7, 0.72), color=INK2, lw=1.0, ms=5)
    arrow(ax, (0.08, 0.6), (0.3, 0.6), color=MUTED, lw=0.9, ms=5)
    ax.text(0.19, 0.65, "$h$", ha="center", va="bottom", fontsize=6.5, color=MUTED)
    ax.text(0.72, 0.52, "easy\naxis", ha="left", va="center", fontsize=5.0, color=MUTED, linespacing=0.9)


def gene(ax, x, y, w, h, label):
    ax.add_patch(mp.FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.035",
                                   fc=FILL, ec=INK2, lw=0.6))
    ax.text(x, y, label, ha="center", va="center", fontsize=6.5, color=INK2)


def draw_toggle(ax, unequal=False):
    gene(ax, 0.27, 0.72, 0.26, 0.16, "A")
    wb, hb = (0.2, 0.12) if unequal else (0.26, 0.16)
    gene(ax, 0.73, 0.72, wb, hb, "B")
    repress(ax, (0.3, 0.81), (0.71, 0.79), -0.45)
    repress(ax, (0.7, 0.65), (0.3, 0.63), -0.45)
    if unequal:
        ax.text(0.73, 0.57, r"$\gamma$", ha="center", va="center", fontsize=6.0, color=AQUA)


def draw_ring(ax):
    c, r = np.array([0.5, 0.72]), 0.19
    pts = [c + r * np.array([np.cos(a), np.sin(a)]) for a in np.radians([90, 0, -90, 180])]
    for i, p in enumerate(pts):
        q = pts[(i + 1) % 4]
        a, b = p + 0.25 * (q - p), p + 0.75 * (q - p)
        repress(ax, tuple(a), tuple(b), -0.25)
    for p in pts:
        ax.add_patch(mp.Circle(tuple(p), 0.055, fc=FILL, ec=INK2, lw=0.6, zorder=3))


def draw_reactor(ax):
    ax.plot([0.2, 0.2, 0.8, 0.8], [0.95, 0.5, 0.5, 0.95], color=INK2, lw=0.8)
    ax.add_patch(mp.Rectangle((0.2, 0.5), 0.6, 0.31, fc=FILL, ec="none"))
    ax.text(0.5, 0.72, r"A + 2X $\rightleftharpoons$ 3X", ha="center", va="center", fontsize=5.2, color=INK2)
    ax.text(0.5, 0.6, r"X $\rightleftharpoons$ B", ha="center", va="center", fontsize=5.2, color=INK2)
    arrow(ax, (0.04, 0.97), (0.28, 0.88), color=MUTED, lw=0.8, ms=4, connectionstyle="arc3,rad=-0.3")
    ax.text(0.07, 0.86, "$b$", ha="center", va="center", fontsize=6.5, color=MUTED)


DRAW = {"pitchfork": draw_landau, "laser": draw_laser, "convection_roll": draw_roll, "stoner_wohlfarth": draw_magnet,
        "toggle": draw_toggle, "repressor_ring4": draw_ring, "schlogl": draw_reactor,
        "toggle_unequal": lambda ax: draw_toggle(ax, unequal=True)}


# the two written states: blue for the state that the bias favours, orange for the other
def state_icon(ax, key, x, y, sign, color, k=1.45):
    s = 0.045 * k
    if key == "pitchfork":
        arrow(ax, (x - sign * s, y), (x + sign * s, y), color=color, lw=1.0, ms=4)
    elif key == "laser":
        t = np.linspace(-1, 1, 40)
        ax.plot(x + s * 1.1 * t, y + sign * 0.03 * k * np.sin(pi * t * 1.5), color=color, lw=0.9)
    elif key == "convection_roll":
        if sign > 0:
            ellipse_arc(ax, (x, y), 0.035 * k, 0.03 * k, 160, -150, color, ms=4)
        else:
            ellipse_arc(ax, (x, y), 0.035 * k, 0.03 * k, 20, 330, color, ms=4)
    elif key == "stoner_wohlfarth":
        arrow(ax, (x, y - sign * 0.035 * k), (x, y + sign * 0.035 * k), color=color, lw=1.0, ms=4)
    elif key in ("toggle", "toggle_unequal"):
        hi, lo = 0.06 * k, 0.02 * k
        a, b = (hi, lo) if sign > 0 else (lo, hi)
        ax.add_patch(mp.Rectangle((x - 0.03 * k, y - 0.03 * k), 0.022 * k, a, color=color, lw=0))
        ax.add_patch(mp.Rectangle((x + 0.008 * k, y - 0.03 * k), 0.022 * k, b, color=color, lw=0))
    elif key == "repressor_ring4":
        for i, (dx, dy) in enumerate(((0, 1), (1, 0), (0, -1), (-1, 0))):
            on = (i % 2 == 0) == (sign > 0)
            ax.add_patch(mp.Circle((x + 0.028 * k * dx, y + 0.028 * k * dy), 0.012 * k, fc=color if on else "white", ec=color,
                                   lw=0.6))
    elif key == "schlogl":
        a, b = 0.025 * k, 0.03 * k
        ax.plot([x - a, x - a, x + a, x + a], [y + 1.2 * b, y - b, y - b, y + 1.2 * b], color=color, lw=0.7)
        ax.add_patch(mp.Rectangle((x - a, y - b), 2 * a, (1.7 if sign > 0 else 0.5) * b, color=color, lw=0))


FIG_W, FIG_H = 180.0, 120.0  # mm
PIC_Y = (0.44, 0.99)  # the vertical range of a pictogram, in units of its width


def mm_axes(fig, x, y, w, h):
    """Axes placed in millimetres from the lower left corner of the figure."""
    return fig.add_axes([x / FIG_W, y / FIG_H, w / FIG_W, h / FIG_H])


def panel_systems(fig, top, left=3.0, width=176.0):
    """A row of tiles: pictogram, name and field, the swept parameter, and the two written states."""
    tw = width / len(MODELS)
    pw = 0.94 * tw
    ph = (PIC_Y[1] - PIC_Y[0]) * pw
    for i, key in enumerate(MODELS):
        m = META[key]
        xl = left + i * tw
        y = top - ph
        pic = mm_axes(fig, xl + 0.03 * tw, y, pw, ph)
        pic.set_xlim(0, 1)
        pic.set_ylim(*PIC_Y)
        pic.axis("off")
        DRAW[key](pic)
        y -= 1.0
        txt = mm_axes(fig, xl, y - 11.5, tw, 11.5)
        txt.set_xlim(0, 1)
        txt.set_ylim(0, 11.5)
        txt.axis("off")
        txt.plot([0.09], [9.0], marker=m["marker"], ms=4.2, color=color_of(key), mec="white", mew=0.4, clip_on=False)
        txt.text(0.16, 9.0, m["name"], ha="left", va="center", fontsize=6.7, fontweight="bold")
        txt.text(0.5, 5.6, m["field"], ha="center", va="center", fontsize=5.8, color=MUTED)
        txt.text(0.5, 2.4, "sweep: " + m["sweep"], ha="center", va="center", fontsize=5.8,
                 color=AQUA if m.get("cusp") else INK2)
        ih = 8.5
        icons = mm_axes(fig, xl, y - 11.5 - ih, tw, ih)
        icons.set_xlim(0, 1)
        icons.set_ylim(0, ih / tw)
        icons.axis("off")
        yc = 0.62 * ih / tw
        for sign, xx, col in ((1, 0.27, BLUE), (-1, 0.73, ORANGE)):
            state_icon(icons, key, xx, yc, sign, col)
            icons.text(xx, 0.13 * ih / tw, m["states"][0 if sign > 0 else 1], ha="center", va="center",
                       fontsize=5.4, color=col)
        icons.text(0.5, yc, "or", ha="center", va="center", fontsize=5.3, color=MUTED)
        if i:
            fig.add_artist(Line2D([xl / FIG_W] * 2, [(y - 11.5 - ih + 0.5) / FIG_H, (top - 0.5) / FIG_H],
                                  color=RULE, lw=0.5))
    return top - ph - 1.0 - 11.5 - 8.5


# ------------------------------------------------------------------------------------------------ data panels
def panel_landscape(ax, data):
    x = np.linspace(-0.66, 0.66, 200)
    ax.plot(x, x ** 4 / 4, color=INK, lw=1.0, zorder=2, label=r"$x^4/4$")
    for key in MODELS:
        d = data["models"][key]
        xs, g = np.asarray(d["x"]), np.asarray(d["g"])
        i0 = int(np.argmin(np.abs(xs)))
        V = -np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(xs))])
        V -= V[i0]
        idx = np.arange(xs.size) - i0
        keep = np.flatnonzero((((np.abs(xs) > 0.25) & (idx % 2 == 0)) | ((np.abs(xs) <= 0.25) & (idx % 6 == 0)))
                              & (np.abs(xs) <= 0.66))
        ax.plot(xs[keep], V[keep], ls="none", marker=META[key]["marker"], ms=3.2, mfc="white",
                mec=color_of(key), mew=0.7, zorder=3)
    ax.set_xlim(-0.75, 0.75)
    ax.set_ylim(-0.003, 0.074)
    ax.set_xlabel(r"critical coordinate, $x = |a_3|^{1/2}\,s$")
    ax.set_ylabel("landscape at the write point, $V(x)$")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.0, 0.012, r"$V = x^4/4$", ha="center", va="bottom", fontsize=6.5)
    # the landscape before, at and after the write point, with the bias that tilts it
    for j, (eps, lab) in enumerate(((-1.0, "before"), (0.0, "at the write point"), (1.0, "after"))):
        ins = ax.inset_axes([0.06 + j * 0.31, 0.7, 0.26, 0.25])
        u = np.linspace(-1.7, 1.7, 100)
        V = -eps * u ** 2 / 2 + u ** 4 / 4 - 0.12 * u
        ins.plot(u, V, color=INK2, lw=0.9)
        um = u[np.argmin(V)]
        ins.add_patch(mp.Circle((um, V.min() + 0.33), 0.2, color=BLUE if eps > 0 else MUTED, zorder=3))
        ins.set_xlim(-1.8, 1.8)
        ins.set_ylim(-0.9, 2.3)
        ins.axis("off")
        ins.text(0, -0.85, lab, ha="center", va="top", fontsize=5.8, color=MUTED)
    ax.annotate("", xy=(0.97, 0.985), xytext=(0.06, 0.985), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.6, mutation_scale=5))
    ax.text(0.5, 1.0, "sweep", transform=ax.transAxes, ha="center", va="bottom", fontsize=5.8, color=MUTED)


def panel_laser(ax, data):
    L = data["laser_sweep"]
    P, E = np.asarray(L["P"]), np.asarray(L["E"])
    Pc = L["P_c"]
    end = E[:, -1]
    order = np.random.default_rng(0).permutation(len(E))  # neither colour drawn over the other
    for e in E[order]:
        ax.plot(P, e, color=BLUE if e[-1] > 0 else ORANGE, lw=0.45, alpha=0.4, zorder=2, rasterized=True)
    pp = np.linspace(Pc, P.max(), 100)
    for s in (1, -1):
        ax.plot(pp, s * np.sqrt(pp - Pc), color=INK, lw=0.8, zorder=4)
    ax.plot([P.min(), Pc], [0, 0], color=INK, lw=0.8, zorder=4)
    ax.plot([Pc, P.max()], [0, 0], color=INK, lw=0.8, ls=(0, (2.5, 2)), zorder=4)
    frac = float(np.mean(end > 0))
    ax.text(Pc + 0.02, -1.42, "threshold", ha="left", va="center", fontsize=5.8, color=MUTED)
    ax.axvline(Pc, color=MUTED, lw=0.5, ls=(0, (1, 2)), zorder=1)
    ax.text(0.45, 1.2, f"{frac:.0%} written as phase 0", ha="left", va="center", fontsize=6.5, color=BLUE)
    ax.text(0.45, -1.2, f"{1 - frac:.0%} as phase π", ha="left", va="center", fontsize=6.5, color=ORANGE)
    ax.text(0.45, 0.98, f"law: {L['target_p']:.0%}", ha="left", va="center", fontsize=5.8, color=MUTED)
    # the bias acts over the linear stage, a Gaussian window of width r^(-1/2) in time around the crossing:
    # +- r^(1/2) in the growth rate of the critical mode, +- r^(1/2) / g in the pump
    w = sqrt(L["rate"]) * abs(np.mean(np.diff(P)) / np.mean(np.diff(L["t"]))) / L["rate"]
    ax.axvspan(Pc - w, Pc + w, color=FILL, lw=0, zorder=0)
    ax.text(Pc, 1.46, "the bias acts here", ha="center", va="center", fontsize=5.8, color=MUTED)
    ax.set_xlim(0.4, P.max())
    ax.set_ylim(-1.52, 1.55)
    ax.set_xlabel("pump $P$, swept upward")
    ax.set_ylabel("field amplitude $E$")
    ax.spines[["top", "right"]].set_visible(False)


def panel_law(ax, data, record):
    th = np.linspace(0, 1.6, 200)
    ax.plot(th, [phi(LAW * t) for t in th], color=INK, lw=1.0, zorder=2)
    for key in MODELS:
        rows = data["models"][key]["law"]
        t = np.array([r["theta"] for r in rows])
        P = np.array([r["measured"] for r in rows])
        e = np.array([r["stderr"] for r in rows])
        ax.errorbar(t, P, yerr=e, ls="none", marker=META[key]["marker"], ms=3.4, mfc="white", mec=color_of(key),
                    mew=0.7, ecolor=color_of(key), elinewidth=0.6, capsize=0, zorder=3)
    ax.set_xlim(0, 1.62)
    ax.set_ylim(0.45, 1.01)
    ax.set_xlabel(r"scaled bias $\theta = h_s / (D_s^{1/2}\, r^{1/4})$")
    ax.set_ylabel("accuracy of the write, $P$")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.06, 0.975, r"$P = \Phi(\pi^{1/4}\,\theta)$", fontsize=7, ha="left", va="center")
    ax.text(0.06, 0.94, "no fitted parameter", fontsize=5.8, ha="left", va="center", color=MUTED)
    # the constant of each system from the law record (25,600 trajectories each)
    ins = ax.inset_axes([0.5, 0.1, 0.47, 0.42])
    for j, key in enumerate(MODELS):
        rec = record["nodes"][key]["symmetric-write"]
        ins.errorbar(rec["constant"], j, xerr=rec["stderr"], marker=META[key]["marker"], ms=3.2, mfc="white",
                     mec=color_of(key), mew=0.7, ecolor=color_of(key), elinewidth=0.6, capsize=0)
    ins.axvline(LAW, color=INK, lw=0.8)
    ins.set_ylim(len(MODELS) - 0.4, -0.6)
    ins.set_xlim(1.26, 1.40)
    ins.set_yticks([])
    ins.set_xticks([1.28, LAW, 1.38], ["1.28", r"$\pi^{1/4}$", "1.38"])
    ins.tick_params(labelsize=5.5, length=2)
    ins.spines[["top", "right", "left"]].set_visible(False)
    ins.set_title("constant of each system", fontsize=5.8, color=MUTED, pad=2)


# ------------------------------------------------------------------------------------------------ data
TARGETS = [0.55, 0.65, 0.75, 0.85, 0.93, 0.98]  # accuracies predicted by the law, one bias each
N_TRAJ = 2000


def target(path):
    """The realization, the derivation and the write point at which derive_symmetric_write tests its law (the
    realization along the sweep through a cusp for an unfolded system), with the signs of the bias it uses."""
    from fieldbridge.memory import codiscovery, spec
    from fieldbridge.memory.discovery import card_rng
    real = spec.load(path)
    got = {}

    def grab(r, nf, rng, **kw):
        got.update(real=r, nf=nf, signs=kw.get("signs", (1,)))
        return {}

    keep = codiscovery.swept_write_check
    codiscovery.swept_write_check = grab
    try:
        row = codiscovery.derive_symmetric_write(real, card_rng(SEED, real.name), check_law=True)
    finally:
        codiscovery.swept_write_check = keep
    return real, row, got


def compute(out=DATA):
    from fieldbridge.memory import analysis as an, construct
    from fieldbridge.memory.discovery import card_rng
    data = {"seed": SEED, "targets": TARGETS, "trajectories": N_TRAJ, "models": {}}
    for name in MODELS:
        t0 = time.time()
        real, row, got = target(ROOT / "examples" / "memory" / f"{name}.json")
        rows = []
        for k, p in enumerate(TARGETS):
            law = construct.swept_write_check(got["real"], got["nf"], card_rng(SEED + 1000 + k, real.name),
                                              target_p=p, n_traj=N_TRAJ, replicates=4, signs=got["signs"])
            rows.append({"target_p": p, "theta": NormalDist().inv_cdf(p) / pi ** 0.25, "measured": law["measured"],
                         "stderr": law["stderr"], "gamma": law["gamma"],
                         "by_sign": [{k_: b[k_] for k_ in ("sign", "measured", "stderr")}
                                     for b in law.get("by_sign", [])]})
        c = row["canonical"]
        data["models"][name] = {"name": real.name, "field": real.spec.get("field"), "word": row["word"],
                                "status": row["status"], "signs": list(got["signs"]), "x": c["x"], "g": c["g"],
                                "write_point": row["write_point"], "law": rows}
        print(name, row["word"], [round(r["measured"], 3) for r in rows], round(time.time() - t0), "s", flush=True)
    # the laser swept through threshold, as in the law test at a predicted accuracy of 0.8
    real, row, got = target(ROOT / "examples" / "memory" / "laser.json")
    nf, r = got["nf"], got["real"]
    rate, gamma, p = 0.02, 0.05, 0.8
    vc, qc = nf["value"], np.asarray(nf["state"], float)
    vec, w = np.asarray(nf["mode"], float), np.asarray(nf["left_mode"], float)
    g, a3 = abs(nf["kappa_slope"]), abs(nf["a3"])
    tgt = qc + 1e3 * an.base_length(r) * vec
    h_per = float(w @ r.write_force(qc, tgt, 1.0))
    t_half = 4.0 / sqrt(rate)
    dv = np.sign(-nf["kappa_slope"]) * rate / g
    path = lambda t: vc - dv * t_half + dv * t
    D = gamma * rate / (a3 * float(w @ w))
    h_s = NormalDist().inv_cdf(p) * sqrt(D * float(w @ w)) * rate ** 0.25 / pi ** 0.25
    q0 = np.repeat(qc[None], 150, axis=0)
    _, ts, vals = an.integrate(r, q0, 2 * t_half, D, np.random.default_rng(SEED), write=(tgt, h_s / h_per),
                               control_path=path, record=lambda q: q[:, 0].copy(), n_record=600)
    data["laser_sweep"] = {"t": ts.tolist(), "P": [float(path(t)) for t in ts],
                           "E": np.round(np.asarray(vals).T, 4).tolist(),
                           "P_c": float(vc), "rate": rate, "gamma": gamma, "target_p": p, "D": float(D),
                           "seed_field": float(h_s / h_per), "mode": vec.tolist()}
    Path(out).write_text(json.dumps(data))
    print("written", out)


def draw(data_path, stem):
    data = json.load(open(data_path))
    record = json.load(open(ROOT / "docs" / "site" / "law_constants.json"))
    fig = plt.figure(figsize=(FIG_W / 25.4, FIG_H / 25.4))
    fig.text(1.0 / FIG_W, 1 - 1.5 / FIG_H, "a", fontsize=9, fontweight="bold", va="top")
    bottom = panel_systems(fig, top=FIG_H - 2.0)
    h = bottom - 22.0
    axes = [mm_axes(fig, x, 11.0, w, h) for x, w in ((14.0, 44.0), (74.0, 44.0), (133.0, 45.0))]
    for ax, letter in zip(axes, "bcd"):
        fig.text((ax.get_position().x0 * FIG_W - 12.0) / FIG_W, (bottom - 3.0) / FIG_H, letter, fontsize=9,
                 fontweight="bold", va="top")
    panel_landscape(axes[0], data)
    if "laser_sweep" in data:
        panel_laser(axes[1], data)
    if all("law" in data["models"][k] for k in MODELS):
        panel_law(axes[2], data, record)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{stem}.{ext}", dpi=300 if ext == "png" else None)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "draw"
    if cmd == "compute":
        compute(sys.argv[2] if len(sys.argv) > 2 else DATA)
    else:
        draw(sys.argv[2] if len(sys.argv) > 2 else DATA, sys.argv[3] if len(sys.argv) > 3 else STEM)
