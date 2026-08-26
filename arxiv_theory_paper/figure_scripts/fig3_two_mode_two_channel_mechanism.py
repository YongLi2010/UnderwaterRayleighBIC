#!/usr/bin/env python3
"""Mechanism-led Fig. 3: a two-mode/two-channel Rayleigh-BIC eigenproblem."""
from pathlib import Path
import csv

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
CSV = (ROOT / "Ni2019_MATLAB" / "results" /
       "fig3_radiation_contributions_180k" / "radiation_contributions.csv")
OUT = ROOT / "arxiv_theory_paper" / "figures"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.4,
    "axes.titlesize": 8.2,
    "axes.titleweight": "semibold",
    "axes.linewidth": 0.65,
    "xtick.major.width": 0.65,
    "ytick.major.width": 0.65,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

NAVY = "#17344F"
TEAL = "#168B8C"
GOLD = "#D68A25"
RED = "#D84B3E"
GREY = "#697278"
LIGHT_GREY = "#E9EDEF"
SUBSTRATE = "#C9D0D4"


def read_csv(path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def label(ax, s, x=-0.07, y=1.045):
    ax.text(x, y, s, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9.5, fontweight="bold", color="#111111")


def clear(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)


def geometry_panel(ax):
    clear(ax); label(ax, "a")
    ax.set_title("Periodic two-cavity eigenproblem", pad=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    y0 = 0.47
    ax.add_patch(Rectangle((0.03, 0.18), 0.94, y0-0.18,
                           facecolor=SUBSTRATE, edgecolor="none"))
    ax.plot([0.03, 0.97], [y0, y0], color=NAVY, lw=1.0)
    centers = [0.15, 0.39, 0.63, 0.87]
    for j, xc in enumerate(centers):
        # identical unit cells, each containing the same wide/narrow pair
        ax.add_patch(Rectangle((xc-0.079, y0-0.24), 0.115, 0.24,
                               facecolor="white", edgecolor=NAVY, lw=0.75))
        ax.add_patch(Rectangle((xc+0.058, y0-0.115), 0.030, 0.115,
                               facecolor="white", edgecolor=NAVY, lw=0.75))
    # Highlight one cell without enclosing the whole drawing in a heavy box.
    ax.plot([0.51, 0.75], [0.135, 0.135], color=RED, lw=1.15)
    ax.plot([0.51, 0.51], [0.13, 0.16], color=RED, lw=1.0)
    ax.plot([0.75, 0.75], [0.13, 0.16], color=RED, lw=1.0)
    ax.text(0.63, 0.085, r"unit cell $a$", ha="center", color=RED)

    # Bloch phase and the two physical exterior channels.
    ax.annotate("", xy=(0.80, 0.77), xytext=(0.34, 0.77),
                arrowprops=dict(arrowstyle="-|>", color=NAVY, lw=1.15))
    ax.text(0.57, 0.82, r"Bloch wavevector $k_x$", ha="center", color=NAVY)
    ax.annotate("", xy=(0.63, 0.91), xytext=(0.63, 0.54),
                arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.0))
    ax.text(0.66, 0.69, r"$n=0$", color=GREY, va="center")
    ax.annotate("", xy=(0.18, 0.57), xytext=(0.48, 0.57),
                arrowprops=dict(arrowstyle="-|>", color=GOLD, lw=1.15))
    ax.text(0.25, 0.61, r"$n=-1$", color=GOLD, ha="center")


def network_panel(ax):
    clear(ax); label(ax, "b")
    ax.set_title(r"Radiation matrix $D(k_x)$", pad=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    nodes_left = [(0.15, 0.68, TEAL, r"$a_{\rm L}$"),
                  (0.15, 0.32, GOLD, r"$a_{\rm S}$")]
    nodes_right = [(0.86, 0.68, GREY, r"$n=0$"),
                   (0.86, 0.32, GOLD, r"$n=-1$")]
    for x, y, color, text in nodes_left:
        ax.add_patch(Circle((x, y), 0.075, facecolor=color,
                            edgecolor="white", lw=0.8))
        ax.text(x, y, text, color="white", ha="center", va="center",
                fontsize=8, fontweight="semibold")
    for x, y, color, text in nodes_right:
        ax.add_patch(Circle((x, y), 0.075, facecolor="white",
                            edgecolor=color, lw=1.35))
        ax.text(x, y, text, color=color, ha="center", va="center",
                fontsize=8, fontweight="semibold")

    # Four couplings are the visual 2x2 matrix; no long equation is needed.
    for yl, color in [(0.68, TEAL), (0.32, GOLD)]:
        for yr in (0.68, 0.32):
            rad = 0.0 if abs(yl-yr)<0.01 else (0.18 if yl>yr else -0.18)
            ax.add_patch(FancyArrowPatch((0.23, yl), (0.78, yr),
                         connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>",
                         mutation_scale=7.5, lw=1.05, color=color, alpha=0.88))
    ax.text(0.15, 0.91, "cavity modes", ha="center", color=NAVY)
    ax.text(0.86, 0.91, "radiation channels", ha="center", color=NAVY)
    ax.text(0.50, 0.08, "four complex couplings", ha="center", color=NAVY)


def phasor_inset(ax, z1, z2):
    ia = inset_axes(ax, width="27%", height="58%", loc="upper right",
                    borderpad=0.7)
    clear(ia); ia.set_aspect("equal")
    m = max(abs(z1), abs(z2))*1.25
    ia.set_xlim(-m, m); ia.set_ylim(-m, m)
    ia.axhline(0, color="#D9DEE1", lw=0.55)
    ia.axvline(0, color="#D9DEE1", lw=0.55)
    for z, color in [(z1, TEAL), (z2, GOLD)]:
        ia.add_patch(FancyArrowPatch((0, 0), (z.real, z.imag),
                     arrowstyle="-|>", mutation_scale=7.5, lw=1.25,
                     color=color, shrinkA=0, shrinkB=0))
    ia.add_patch(Circle((0, 0), m*0.045, facecolor=RED,
                        edgecolor="white", lw=0.35, zorder=5))


def channel_panel(fig, spec, data, prefix, title, panel_letter, k_bic):
    outer = fig.add_subplot(spec)
    clear(outer); label(outer, panel_letter, x=-0.075, y=1.045)
    outer.set_title(title, pad=5)
    sub = spec.subgridspec(2, 1, hspace=0.10)
    ar = fig.add_subplot(sub[0, 0])
    ai = fig.add_subplot(sub[1, 0], sharex=ar)
    k = data["kappa"]
    l = data[f"{prefix}_large_real"] + 1j*data[f"{prefix}_large_imag"]
    s = data[f"{prefix}_small_real"] + 1j*data[f"{prefix}_small_imag"]
    ib = int(np.argmin(abs(k-k_bic)))

    for ax, lv, sv, ylabel in [(ar, l.real, s.real, "Re"),
                               (ai, l.imag, s.imag, "Im")]:
        ax.plot(k, lv, color=TEAL, lw=1.45)
        ax.plot(k, sv, color=GOLD, lw=1.45)
        ax.axvline(k_bic, color=RED, lw=0.8, ls=(0, (2.5, 2)))
        ax.axhline(0, color="#C8CED1", lw=0.55)
        ax.scatter([k_bic, k_bic], [lv[ib], sv[ib]], s=17,
                   c=[TEAL, GOLD], edgecolors="white", linewidths=0.45,
                   zorder=5)
        lim = max(abs(np.r_[lv, sv]))*1.14
        ax.set_ylim(-lim, lim)
        ax.set_ylabel(ylabel, rotation=0, rotation_mode="anchor",
                      labelpad=10, color=NAVY)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=6.3, pad=1.5)
    ar.tick_params(labelbottom=False)
    ai.set_xlabel(r"Bloch wavevector $k_xa/2\pi$")
    ar.text(k_bic, ar.get_ylim()[1]*0.88, "BIC", ha="center", va="top",
            color=RED, fontsize=6.5)
    # The inset encodes the invariant statement: the two complex channel
    # contributions are diametrically opposite at the BIC.
    phasor_inset(ar, l[ib], s[ib])
    return ar, ai


def main():
    data = read_csv(CSV)
    k_bic = data["kappa"][np.argmin(data["b0_sum_abs"] + data["bm1_sum_abs"])]

    fig = plt.figure(figsize=(7.0, 5.0), facecolor="white")
    gs = fig.add_gridspec(2, 2, left=0.065, right=0.985, bottom=0.085,
                          top=0.95, wspace=0.25, hspace=0.32,
                          height_ratios=[0.86, 1.14])
    axa = fig.add_subplot(gs[0, 0]); geometry_panel(axa)
    axb = fig.add_subplot(gs[0, 1]); network_panel(axb)
    channel_panel(fig, gs[1, 0], data, "b0", r"$n=0$ channel",
                  "c", k_bic)
    channel_panel(fig, gs[1, 1], data, "bm1", r"$n=-1$ Rayleigh channel",
                  "d", k_bic)

    OUT.mkdir(parents=True, exist_ok=True)
    stem = OUT / "fig3_two_mode_two_channel_mechanism"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(stem)


if __name__ == "__main__":
    main()
