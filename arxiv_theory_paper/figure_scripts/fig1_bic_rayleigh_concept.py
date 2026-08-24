"""PRL Fig. 1: conventional BIC, Rayleigh anomaly, and strict Rayleigh BIC."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, Polygon, Rectangle
import numpy as np


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 8.2,
        "axes.linewidth": 0.8,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "legend.frameon": False,
        "mathtext.fontset": "dejavusans",
    }
)

BLUE = "#2166AC"
BLUE_LIGHT = "#DCEAF5"
ORANGE = "#D97706"
ORANGE_LIGHT = "#FCE8CB"
RED = "#C7352B"
GREEN = "#287A4D"
GREEN_LIGHT = "#DDEFE5"
INK = "#20242A"
MID = "#6B7280"
PALE = "#F5F6F7"

# Shared layout grid: identical anchors across all three conceptual panels.
TITLE_Y = 0.955
KEY_Y = 0.750
SURFACE_Y = 0.220
OPEN_CANCEL_Y = 0.660


def setup_panel(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")


def panel_label(ax, letter):
    ax.text(
        0.01,
        0.985,
        f"({letter})",
        ha="left",
        va="top",
        fontsize=10.5,
        fontweight="bold",
        color=INK,
    )


def title(ax, text):
    ax.text(0.50, TITLE_Y, text, ha="center", va="top", fontsize=9.0, fontweight="bold", color=INK)


def arrow(ax, xy0, xy1, color=INK, lw=1.4, style="-|>", mutation=9, connection="arc3"):
    patch = FancyArrowPatch(
        xy0,
        xy1,
        arrowstyle=style,
        mutation_scale=mutation,
        lw=lw,
        color=color,
        connectionstyle=connection,
        shrinkA=0,
        shrinkB=0,
    )
    ax.add_patch(patch)
    return patch


def cross(ax, x, y, size=0.032, color=RED, lw=1.7):
    ax.plot([x - size, x + size], [y - size, y + size], color=color, lw=lw, solid_capstyle="round")
    ax.plot([x - size, x + size], [y + size, y - size], color=color, lw=lw, solid_capstyle="round")


def periodic_boundary(ax, y0=0.27, cavity=True):
    bottom = 0.08
    ax.add_patch(Rectangle((0.05, bottom), 0.90, y0 - bottom, facecolor="#E7EAED", edgecolor="none"))
    ax.plot([0.05, 0.95], [y0, y0], color=INK, lw=1.6)
    for x in [0.10, 0.90]:
        ax.plot([x, x], [bottom, y0 + 0.075], color=MID, lw=0.7, ls=(0, (2.0, 2.0)))
    # Split the dimension line around the label so no stroke crosses the glyph.
    arrow(ax, (0.445, 0.045), (0.10, 0.045), MID, 0.8, style="-|>", mutation=6)
    arrow(ax, (0.555, 0.045), (0.90, 0.045), MID, 0.8, style="-|>", mutation=6)
    ax.text(0.50, 0.045, r"$a$", ha="center", va="center", fontsize=7.5, color=INK)
    if cavity:
        ax.add_patch(Rectangle((0.34, bottom), 0.16, y0 - bottom, facecolor=BLUE_LIGHT, edgecolor=BLUE, lw=1.1))
        ax.add_patch(Rectangle((0.62, bottom), 0.09, y0 - bottom, facecolor=ORANGE_LIGHT, edgecolor=ORANGE, lw=1.1))


def draw_conventional_bic(ax):
    setup_panel(ax)
    panel_label(ax, "a")
    title(ax, "Conventional BIC")

    ax.plot([0.10, 0.90], [SURFACE_Y, SURFACE_Y], color=INK, lw=1.4)
    ax.add_patch(Ellipse((0.50, 0.15), 0.40, 0.18, facecolor="#EEF3F7", edgecolor=INK, lw=1.0))
    ax.add_patch(Circle((0.41, 0.16), 0.045, facecolor=BLUE, edgecolor="white", lw=0.8))
    ax.add_patch(Circle((0.59, 0.16), 0.045, facecolor=ORANGE, edgecolor="white", lw=0.8))
    ax.text(0.41, 0.16, r"$q_1$", color="white", ha="center", va="center", fontsize=7.5)
    ax.text(0.59, 0.16, r"$q_2$", color="white", ha="center", va="center", fontsize=7.5)

    arrow(ax, (0.41, 0.22), (0.47, OPEN_CANCEL_Y - 0.01), BLUE, 1.7, connection="arc3,rad=-0.18")
    arrow(ax, (0.59, 0.22), (0.53, OPEN_CANCEL_Y - 0.01), ORANGE, 1.7, connection="arc3,rad=0.18")
    cross(ax, 0.50, OPEN_CANCEL_Y, size=0.030)
    ax.text(0.50, KEY_Y, r"$A_{\mathrm{open}}=0$", ha="center", va="center", fontsize=9.0, color=INK)


def draw_rayleigh_anomaly(ax):
    setup_panel(ax)
    panel_label(ax, "b")
    title(ax, "Rayleigh anomaly")
    periodic_boundary(ax, y0=SURFACE_Y, cavity=False)

    # One generic groove per dashed unit cell; the neighboring cells are implied
    # by the period boundaries and the dimension a.
    ax.add_patch(Rectangle((0.44, 0.08), 0.12, SURFACE_Y - 0.08,
                           facecolor="#F8F9FA", edgecolor="none"))
    ax.plot([0.44, 0.44, 0.56, 0.56], [SURFACE_Y, 0.10, 0.10, SURFACE_Y],
            color=MID, lw=1.0)

    for y in [0.38, 0.50, 0.62]:
        ax.plot([0.12, 0.88], [y, y], color=ORANGE, lw=0.9, alpha=0.22)
    arrow(ax, (0.18, 0.50), (0.84, 0.50), ORANGE, 2.0, mutation=10)
    ax.text(0.50, KEY_Y, r"$k_{y,\mathrm{RA}}=0$", ha="center", va="center", fontsize=9.0, color=INK)
    ax.text(0.50, 0.31, r"$J_y=0$", ha="center", va="center", fontsize=8.8, color=INK)
    ax.text(0.50, 0.44, r"$A_{\mathrm{RA}}\ne0$", ha="center", va="center", fontsize=8.8, color=RED)


def draw_rayleigh_bic(ax):
    setup_panel(ax)
    panel_label(ax, "c")
    title(ax, "Rayleigh BIC")
    periodic_boundary(ax, y0=SURFACE_Y, cavity=True)

    ax.add_patch(Polygon([[0.34, 0.22], [0.50, 0.22], [0.47, 0.34], [0.37, 0.34]], closed=True,
                         facecolor=BLUE, edgecolor="none", alpha=0.88))
    ax.add_patch(Polygon([[0.62, 0.22], [0.71, 0.22], [0.69, 0.31], [0.64, 0.31]], closed=True,
                         facecolor=ORANGE, edgecolor="none", alpha=0.88))
    ax.text(0.42, 0.145, r"$q_1$", ha="center", va="center", color=BLUE, fontsize=7.5, fontweight="bold")
    ax.text(0.665, 0.145, r"$q_2$", ha="center", va="center", color=ORANGE, fontsize=7.5, fontweight="bold")

    arrow(ax, (0.43, 0.34), (0.57, OPEN_CANCEL_Y - 0.01), BLUE, 1.7, connection="arc3,rad=-0.12")
    arrow(ax, (0.665, 0.31), (0.61, OPEN_CANCEL_Y - 0.01), ORANGE, 1.7, connection="arc3,rad=0.14")
    cross(ax, 0.59, OPEN_CANCEL_Y, size=0.028)
    ax.text(0.59, KEY_Y, r"$A_{\mathrm{open}}=0$", ha="center", va="center", fontsize=8.8, color=INK)

    arrow(ax, (0.43, 0.32), (0.83, 0.47), BLUE, 1.7, connection="arc3,rad=0.13")
    arrow(ax, (0.665, 0.30), (0.83, 0.47), ORANGE, 1.7, connection="arc3,rad=-0.13")
    cross(ax, 0.85, 0.48, size=0.027)
    ax.text(0.78, 0.57, r"$A_{\mathrm{RA}}=0$", ha="center", va="center", fontsize=8.8, color=INK)


def main():
    fig = plt.figure(figsize=(7.10, 2.35), facecolor="white")
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.0], wspace=0.045)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    draw_conventional_bic(axes[0])
    draw_rayleigh_anomaly(axes[1])
    draw_rayleigh_bic(axes[2])

    fig.subplots_adjust(left=0.012, right=0.988, bottom=0.035, top=0.985)
    fig.canvas.draw()
    positions = [ax.get_position() for ax in axes]
    assert max(p.width for p in positions) - min(p.width for p in positions) < 1e-10
    assert max(p.height for p in positions) - min(p.height for p in positions) < 1e-10
    assert max(p.y0 for p in positions) - min(p.y0 for p in positions) < 1e-10
    out_dir = Path(__file__).resolve().parents[1] / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / "fig1_bic_rayleigh_concept"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.02)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.02)
    fig.savefig(stem.with_suffix(".png"), dpi=400, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


if __name__ == "__main__":
    main()
