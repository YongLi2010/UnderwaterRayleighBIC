"""Conceptual PRL/Nature-style Fig. 1: BIC, Rayleigh anomaly, Rayleigh BIC.

The artwork is deliberately schematic rather than data-like. All objects are
vector artists and remain editable in the exported PDF/SVG.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, Ellipse, FancyArrowPatch
import numpy as np


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8.5,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "mathtext.fontset": "dejavusans",
        "lines.solid_capstyle": "round",
        "lines.solid_joinstyle": "round",
    }
)

# Restrained, color-blind-safe palette. Red is reserved for the radiation zero.
INK = "#20252B"
NAVY = "#173A63"
TEAL = "#248796"
ORANGE = "#D98127"
RED = "#C43B32"
GRAY = "#AAB2B9"
LIGHT_GRAY = "#E7EAED"
PALE_BLUE = "#E7EFF6"


def setup(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")


def header(ax, label, title):
    ax.text(0.025, 0.965, f"({label})", ha="left", va="top", fontsize=10.5,
            fontweight="bold", color=INK)
    ax.text(0.50, 0.925, title, ha="center", va="top", fontsize=9.4,
            fontweight="bold", color=INK)


def arrow(ax, start, end, color, lw=1.35, scale=8, connection="arc3"):
    patch = FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=scale, lw=lw,
        color=color, connectionstyle=connection, shrinkA=0, shrinkB=0,
    )
    ax.add_patch(patch)
    return patch


def zero_node(ax, x, y, radius=0.019):
    ax.add_patch(Circle((x, y), radius, facecolor="white", edgecolor=RED, lw=1.55, zorder=9))
    ax.add_patch(Circle((x, y), radius * 0.22, facecolor=RED, edgecolor="none", zorder=10))


def localized_state(ax, center=(0.50, 0.47), width=0.44, height=0.34):
    """Nested, phase-alternating contours suggesting a localized eigenfield."""
    cx, cy = center
    for scale, color, alpha in [
        (1.00, PALE_BLUE, 0.78),
        (0.77, "#BFD3E4", 0.52),
        (0.55, TEAL, 0.24),
        (0.34, NAVY, 0.82),
    ]:
        ax.add_patch(Ellipse((cx, cy), width * scale, height * scale,
                             facecolor=color, edgecolor="none", alpha=alpha, zorder=2))
    ax.plot([cx, cx], [cy - 0.105, cy + 0.105], color="white", lw=0.75, alpha=0.9, zorder=8)
    ax.add_patch(Ellipse((cx - 0.065, cy), 0.095, 0.13, facecolor=NAVY,
                         edgecolor="white", lw=0.55, alpha=0.92, zorder=7))
    ax.add_patch(Ellipse((cx + 0.065, cy), 0.095, 0.13, facecolor=TEAL,
                         edgecolor="white", lw=0.55, alpha=0.92, zorder=7))


def periodic_interface(ax, y=0.225):
    """Abstract periodic interface, intentionally not a cavity geometry."""
    ax.plot([0.055, 0.945], [y, y], color=INK, lw=1.35, zorder=8)
    xs = np.linspace(0.095, 0.905, 9)
    for x in xs:
        ax.add_patch(Circle((x, y), 0.018, facecolor="white", edgecolor=NAVY, lw=1.05, zorder=9))
        ax.add_patch(Circle((x, y), 0.0062, facecolor=NAVY, edgecolor="none", zorder=10))


def conventional_bic(ax):
    setup(ax)
    header(ax, "a", "Conventional BIC")

    # Faint continuum wavefronts: the state lies inside them but does not feed them.
    for w, h, alpha in [(0.69, 0.55, 0.60), (0.84, 0.68, 0.42)]:
        ax.add_patch(Arc((0.50, 0.47), w, h, theta1=20, theta2=160,
                         color=GRAY, lw=0.8, alpha=alpha, zorder=0))
        ax.add_patch(Arc((0.50, 0.47), w, h, theta1=200, theta2=340,
                         color=GRAY, lw=0.8, alpha=alpha, zorder=0))
    localized_state(ax, center=(0.50, 0.48), width=0.45, height=0.34)

    # Two radiation pathways reach the same continuum port in opposite phase.
    arrow(ax, (0.425, 0.52), (0.49, 0.755), TEAL, connection="arc3,rad=-0.22")
    arrow(ax, (0.575, 0.52), (0.51, 0.755), ORANGE, connection="arc3,rad=0.22")
    zero_node(ax, 0.50, 0.765)
    ax.text(0.50, 0.105, r"$A_{\rm open}=0$", ha="center", va="center",
            fontsize=9.1, color=INK)


def rayleigh_anomaly(ax):
    setup(ax)
    header(ax, "b", "Rayleigh anomaly")
    periodic_interface(ax)

    # ky=0 field: vertical phase fronts persist through the exterior instead of decaying.
    xvals = np.linspace(0.10, 0.90, 7)
    for j, x in enumerate(xvals):
        color = TEAL if j % 2 == 0 else NAVY
        alpha = 0.13 if j % 2 == 0 else 0.075
        ax.add_patch(Ellipse((x, 0.515), 0.105, 0.55, facecolor=color,
                             edgecolor="none", alpha=alpha, zorder=0))
        ax.plot([x, x], [0.275, 0.755], color=color, lw=0.75, alpha=0.40, zorder=1)
    arrow(ax, (0.16, 0.535), (0.86, 0.535), TEAL, lw=1.75, scale=9)
    ax.text(0.50, 0.805, r"$k_{y,\rm RA}=0$", ha="center", va="center",
            fontsize=9.0, color=INK)
    ax.text(0.50, 0.105, r"$J_y=0,\quad A_{\rm RA}\ne0$", ha="center", va="center",
            fontsize=9.0, color=INK)


def rayleigh_bic(ax):
    setup(ax)
    header(ax, "c", "Rayleigh BIC")
    periodic_interface(ax)

    # A surface-localized eigenfield replaces the non-decaying grazing field.
    localized_state(ax, center=(0.50, 0.405), width=0.46, height=0.31)

    # Both independent exterior amplitudes terminate at zeros.
    arrow(ax, (0.58, 0.42), (0.825, 0.53), TEAL, lw=1.35,
          connection="arc3,rad=-0.12")
    zero_node(ax, 0.845, 0.54)
    arrow(ax, (0.47, 0.49), (0.50, 0.735), ORANGE, lw=1.35,
          connection="arc3,rad=0.08")
    zero_node(ax, 0.505, 0.755)

    ax.text(0.50, 0.105, r"$A_{\rm RA}=A_{\rm open}=0$", ha="center", va="center",
            fontsize=9.0, color=INK)


def main():
    # 183 mm is the full PRL two-column width; the short aspect ratio keeps it conceptual.
    fig = plt.figure(figsize=(7.2047, 2.5984), facecolor="white")
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1], wspace=0.035)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    conventional_bic(axes[0])
    rayleigh_anomaly(axes[1])
    rayleigh_bic(axes[2])

    fig.canvas.draw()
    for i in [0, 1]:
        p_left, p_right = axes[i].get_position(), axes[i + 1].get_position()
        x = (p_left.x1 + p_right.x0) / 2
        fig.add_artist(plt.Line2D([x, x], [0.115, 0.90], transform=fig.transFigure,
                                  color=LIGHT_GRAY, lw=0.65))

    fig.subplots_adjust(left=0.012, right=0.988, bottom=0.025, top=0.985)
    out = Path(__file__).resolve().parents[1] / "figures"
    out.mkdir(parents=True, exist_ok=True)
    stem = out / "fig1_bic_rayleigh_concept_nature"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.015)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.015)
    fig.savefig(stem.with_suffix(".png"), dpi=450, bbox_inches="tight", pad_inches=0.015)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pad_inches=0.015)
    plt.close(fig)


if __name__ == "__main__":
    main()
