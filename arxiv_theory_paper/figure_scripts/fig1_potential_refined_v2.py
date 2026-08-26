"""Conceptual Fig. 1 refined from the author's original visual grammar.

Panel b deliberately uses a one-site periodic lattice: a Rayleigh anomaly needs
periodicity, not a two-resonator basis.  Panel c introduces two sites per cell
only where the two independent radiation amplitudes must be cancelled.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, FancyArrowPatch, FancyBboxPatch
import numpy as np


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 8.2,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    # Keep prose, Greek symbols, numerals, and mathematical annotations in one
    # Arial family so that off-Γ titles do not switch visual grammars.
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.sf": "Arial",
    "mathtext.cal": "Arial",
    "mathtext.tt": "Arial",
    "lines.solid_capstyle": "round",
    "lines.solid_joinstyle": "round",
})

INK = "#20252B"
NAVY = "#17528A"
BLUE = "#236FA0"
ORANGE = "#E0781B"
RED = "#D13B32"
GRAY = "#7F8890"
LIGHT = "#D6DCE0"
PALE = "#F5F7F8"
HERO_PALE = "#EEF4F8"
HERO_EDGE = "#9CB8CC"
SITE_R = 0.026

def setup(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")


def arrow(ax, start, end, color=INK, lw=1.25, scale=8, alpha=1, zorder=4):
    p = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=scale,
                        linewidth=lw, color=color, alpha=alpha,
                        shrinkA=0, shrinkB=0, zorder=zorder)
    ax.add_patch(p)
    return p


def double_arrow(ax, x1, x2, y, color=NAVY):
    arrow(ax, (x1, y), (x2, y), color, 0.75, 6.0)
    arrow(ax, (x2, y), (x1, y), color, 0.75, 6.0)


def sphere(ax, x, y, radius, color):
    """Vector-only shaded sphere, close to the original visual language."""
    ax.add_patch(Circle((x, y), radius, facecolor=color, edgecolor="white",
                        linewidth=0.7, zorder=5))
    for frac, alpha in [(0.72, 0.12), (0.47, 0.12), (0.25, 0.16)]:
        ax.add_patch(Circle((x - radius * 0.25, y + radius * 0.28),
                            radius * frac, facecolor="white", edgecolor="none",
                            alpha=alpha, zorder=6))


def wave_trace(ax, start, end, color, alpha=1.0, cycles=8, amplitude=0.006):
    """Draw one clean propagation vector without a superposed wavy trajectory."""
    arrow(ax, start, end, color, lw=1.15, scale=7.7, alpha=alpha, zorder=4)


def red_cross(ax, xy, size=0.022):
    x, y = xy
    ax.plot([x - size, x + size], [y - size, y + size], color=RED, lw=2.0, zorder=12)
    ax.plot([x - size, x + size], [y + size, y - size], color=RED, lw=2.0, zorder=12)


def title(ax, letter, name, subtitle=None, hero=False):
    ax.text(0.02, 0.975, letter, ha="left", va="top", fontsize=9.8,
            fontweight="bold", color=INK)
    ax.text(0.50, 0.97, name, ha="center", va="top", fontsize=9.2,
            fontweight="bold", color=NAVY if hero else INK)
    if subtitle:
        ax.text(0.50, 0.905, subtitle, ha="center", va="top",
                fontsize=7.3, color=GRAY)


def criterion(ax, formula, hero=False):
    box = FancyBboxPatch((0.215, 0.785), 0.57, 0.072,
                         boxstyle="round,pad=0.010,rounding_size=0.016",
                         facecolor=HERO_PALE if hero else PALE,
                         edgecolor=HERO_EDGE if hero else LIGHT,
                         linewidth=0.95 if hero else 0.8, zorder=0)
    ax.add_patch(box)
    ax.text(0.50, 0.821, formula, ha="center", va="center",
            fontsize=8.7, color=INK)


def one_site_lattice(ax, y=0.150):
    ax.plot([0.075, 0.925], [y, y], color=GRAY, lw=0.75)
    # Same lattice and same blue site as panel c; panel c only adds a second basis site.
    sites = [0.31, 0.69]
    for x in sites:
        sphere(ax, x, y + SITE_R, SITE_R, NAVY)
    ax.text(0.035, y + SITE_R, "…", ha="center", va="center", fontsize=10.0, color=GRAY)
    ax.text(0.965, y + SITE_R, "…", ha="center", va="center", fontsize=10.0, color=GRAY)
    double_arrow(ax, sites[0], sites[1], 0.078)
    ax.text((sites[0] + sites[1]) / 2, 0.045, r"$a$", ha="center",
            va="center", fontsize=7.5, color=INK)


def two_site_lattice(ax, y=0.150):
    ax.plot([0.075, 0.925], [y, y], color=GRAY, lw=0.75)
    blue = [0.31, 0.69]
    # The second basis site is inserted between equivalent blue sites;
    # neither the period nor the original lattice is changed.
    orange = [0.12, 0.50, 0.88]
    for x in blue:
        sphere(ax, x, y + SITE_R, SITE_R, NAVY)
    for x in orange:
        sphere(ax, x, y + SITE_R, SITE_R, ORANGE)
    ax.text(0.035, y + SITE_R, "…", ha="center", va="center", fontsize=10.0, color=GRAY)
    ax.text(0.965, y + SITE_R, "…", ha="center", va="center", fontsize=10.0, color=GRAY)
    double_arrow(ax, blue[0], blue[1], 0.078)
    ax.text((blue[0] + blue[1]) / 2, 0.045, r"$a$", ha="center",
            va="center", fontsize=7.5, color=INK)


def channel_diagram(ax, cancelled=False):
    """Minimal off-Gamma topology: a slanted n=0 order and grazing n=-1 order."""
    # Compact the vertical story: criterion -> channel topology -> lattice.
    ox, oy = 0.50, 0.370
    rx, ry = 0.285, 0.270
    axis_overhang = 0.045
    x_axis_left = ox - rx - axis_overhang
    x_axis_right = ox + rx + axis_overhang
    y_axis_top = oy + ry + axis_overhang

    # The pale semicircle is the real-radiation boundary |kx,n|=k.
    ax.add_patch(Arc((ox, oy), 2 * rx, 2 * ry, theta1=0, theta2=180,
                     color=LIGHT, lw=0.9, zorder=0))
    # Single continuous arrow objects prevent seams and make the x/y axis
    # overhang and arrowheads optically identical.
    arrow(ax, (x_axis_left, oy), (x_axis_right, oy), INK, 0.85, 6.2)
    arrow(ax, (ox, oy), (ox, y_axis_top), INK, 0.85, 6.2)
    ax.text(x_axis_right + 0.004, oy - 0.036, r"$k_x$",
            ha="center", va="top", fontsize=7.5)
    ax.text(ox - 0.024, y_axis_top + 0.005, r"$k_y$", ha="right", va="center", fontsize=7.5)

    # A small generic off-Gamma tilt is used conceptually; no design angle is encoded.
    u0 = 0.16
    v0 = np.sqrt(1 - u0 * u0)
    p0 = (ox + rx * u0, oy + ry * v0)
    pm = (ox - rx, oy)

    alpha = 0.22 if cancelled else 1.0
    wave_trace(ax, (ox, oy), p0, BLUE, alpha=alpha, amplitude=0.0055)
    wave_trace(ax, (ox, oy), pm, ORANGE, alpha=alpha, amplitude=0.0055)
    ax.add_patch(Circle((ox, oy), 0.009, facecolor=INK, edgecolor="none", zorder=9))

    ax.text(p0[0] + 0.024, p0[1] + 0.012, r"$n=0$", color=BLUE,
            fontsize=7.6, ha="left", va="bottom")
    # Place the grazing-order label below its arrow, away from the radiation arc.
    ax.text(ox - 0.57 * rx, oy - 0.047, r"$n=-1$", color=ORANGE,
            fontsize=7.6, ha="center", va="top")
    if cancelled:
        red_cross(ax, (ox + 0.58 * (p0[0] - ox), oy + 0.58 * (p0[1] - oy)), 0.020)
        red_cross(ax, (ox + 0.58 * (pm[0] - ox), oy), 0.020)


def panel_a(ax):
    setup(ax)
    title(ax, "a", "Conventional BIC")
    criterion(ax, r"$A_{\rm rad}=0$")

    # The grey substrate lies below the conceptual scattering degrees of
    # freedom; each site is tangent to the upper surface.
    source_y = 0.160 + SITE_R
    sphere(ax, 0.27, source_y, SITE_R, NAVY)
    sphere(ax, 0.73, source_y, SITE_R, ORANGE)
    ax.plot([0.12, 0.88], [0.160, 0.160], color=GRAY, lw=0.75)
    merge = (0.50, 0.495)
    wave_trace(ax, (0.27, source_y + SITE_R), merge, BLUE, cycles=7)
    wave_trace(ax, (0.73, source_y + SITE_R), merge, ORANGE, cycles=7)
    ax.add_patch(Circle(merge, 0.008, facecolor=INK, edgecolor="none", zorder=9))

    # Suppressed radiation channel above the cancellation point.
    arrow(ax, merge, (0.50, 0.735), GRAY, 1.05, 7.2)
    red_cross(ax, (0.50, 0.610), 0.021)


def panel_b(ax):
    setup(ax)
    title(ax, "b", "Rayleigh anomaly  off-Γ")
    criterion(ax, r"$A_0\ne0,\quad A_{-1}\ne0$")
    channel_diagram(ax, cancelled=False)
    one_site_lattice(ax)


def panel_c(ax):
    setup(ax)
    title(ax, "c", "Rayleigh BIC  off-Γ", hero=True)
    criterion(ax, r"$A_0=0,\quad A_{-1}=0$", hero=True)
    channel_diagram(ax, cancelled=True)
    two_site_lattice(ax)


def main():
    fig = plt.figure(figsize=(7.2047, 3.42), facecolor="white")
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1], wspace=0.012)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    panel_a(axes[0])
    panel_b(axes[1])
    panel_c(axes[2])
    fig.subplots_adjust(left=0.012, right=0.992, bottom=0.025, top=0.985)

    out = Path(__file__).resolve().parents[1] / "figures"
    out.mkdir(parents=True, exist_ok=True)
    stem = out / "fig1_potential_refined_v2"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.015)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.015)
    fig.savefig(stem.with_suffix(".png"), dpi=450, bbox_inches="tight", pad_inches=0.015)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pad_inches=0.015)
    plt.close(fig)


if __name__ == "__main__":
    main()
