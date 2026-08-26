"""Corrected conceptual Fig. 1 for the 7.10-degree Rayleigh-BIC design.

Figure contract
---------------
Core conclusion: a Rayleigh anomaly and a Rayleigh BIC share the same Floquet
channel topology, but only the BIC cancels both the finite-flux n=0 amplitude
and the grazing n=-1 pressure amplitude.
Archetype: schematic-led conceptual triptych.
Evidence roles: (a) conventional radiation cancellation; (b) channel opening;
(c) simultaneous channel zeros.
Geometry: kappa=0.1100001053, Omega=1-kappa, theta=7.09966 degrees.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, Ellipse, FancyArrowPatch
import numpy as np


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 8.2,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "mathtext.fontset": "dejavusans",
        "lines.solid_capstyle": "round",
        "lines.solid_joinstyle": "round",
    }
)

INK = "#22272D"
NAVY = "#1E4D7A"
BLUE = "#2C78A0"
TEAL = "#2B8C8E"
ORANGE = "#D97A24"
RED = "#C73B32"
MID = "#7E878F"
LIGHT = "#D9DEE2"
PALE = "#EDF2F5"
PALE_BLUE = "#DCE9F2"

# Verified final root used only to set the physically correct channel angles.
KAPPA = 0.110000105311808
OMEGA = 1.0 - KAPPA
THETA_DEG = np.degrees(np.arcsin(KAPPA / OMEGA))


def setup(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")


def header(ax, letter, title, subtitle=None):
    ax.text(0.018, 0.982, letter, ha="left", va="top", fontsize=10.0,
            fontweight="bold", color=INK)
    ax.text(0.50, 0.974, title, ha="center", va="top", fontsize=9.3,
            fontweight="bold", color=INK)
    if subtitle:
        ax.text(0.50, 0.903, subtitle, ha="center", va="top", fontsize=7.4,
                color=MID)


def arrow(ax, start, end, color, lw=1.45, scale=8.5, alpha=1.0, zorder=5):
    p = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=scale,
                        linewidth=lw, color=color, alpha=alpha,
                        shrinkA=0, shrinkB=0, zorder=zorder)
    ax.add_patch(p)
    return p


def zero_mark(ax, xy, radius=0.018):
    ax.add_patch(Circle(xy, radius, facecolor="white", edgecolor=RED,
                        linewidth=1.45, zorder=12))
    ax.plot([xy[0] - 0.010, xy[0] + 0.010],
            [xy[1] - 0.010, xy[1] + 0.010],
            color=RED, lw=1.15, zorder=13)


def wave_arrow(ax, start, end, color, alpha=1.0, zorder=6):
    """Arrow plus a sparse sinusoidal phase trace."""
    x0, y0 = start
    x1, y1 = end
    vec = np.array([x1 - x0, y1 - y0], dtype=float)
    length = np.hypot(*vec)
    tangent = vec / length
    normal = np.array([-tangent[1], tangent[0]])
    s = np.linspace(0.11, 0.78, 90)
    base = np.array([x0, y0])[:, None] + vec[:, None] * s
    wiggle = normal[:, None] * (0.006 * np.sin(2 * np.pi * 8 * s))[None, :]
    path = base + wiggle
    ax.plot(path[0], path[1], color=color, lw=0.85, alpha=alpha, zorder=zorder)
    arrow(ax, start, end, color, lw=1.35, scale=8.2, alpha=alpha, zorder=zorder + 1)


def localized_mode(ax, center=(0.50, 0.49)):
    cx, cy = center
    for w, h, color, alpha in [
        (0.48, 0.34, PALE, 0.95),
        (0.37, 0.27, PALE_BLUE, 0.90),
        (0.26, 0.19, "#A8C9D8", 0.68),
    ]:
        ax.add_patch(Ellipse((cx, cy), w, h, facecolor=color,
                             edgecolor="none", alpha=alpha, zorder=1))
    ax.add_patch(Ellipse((cx - 0.067, cy), 0.105, 0.145, facecolor=NAVY,
                         edgecolor="white", linewidth=0.6, zorder=3))
    ax.add_patch(Ellipse((cx + 0.067, cy), 0.105, 0.145, facecolor=TEAL,
                         edgecolor="white", linewidth=0.6, zorder=3))
    ax.plot([cx, cx], [cy - 0.105, cy + 0.105], color="white", lw=0.75, zorder=4)


def periodic_dimer(ax, y=0.122):
    """Two complete unit cells; boundaries enclose one blue-orange basis."""
    ax.plot([0.045, 0.955], [y, y], color=MID, lw=0.75, zorder=1)
    # Equivalent blue sites are separated by exactly one period.
    blue_x = [0.17, 0.63]
    orange_x = [0.36, 0.82]
    for x in blue_x:
        ax.add_patch(Circle((x, y + 0.018), 0.030, facecolor=NAVY,
                            edgecolor="white", linewidth=0.65, zorder=5))
    for x in orange_x:
        ax.add_patch(Circle((x, y + 0.018), 0.024, facecolor=ORANGE,
                            edgecolor="white", linewidth=0.65, zorder=5))

    # Boundaries lie between neighbouring dimers, not between the two sites.
    for xb in [0.035, 0.495, 0.955]:
        ax.plot([xb, xb], [0.075, 0.205], color=LIGHT, lw=0.75,
                linestyle=(0, (3, 2)), zorder=0)

    # Period is measured between equivalent blue sites.
    ydim = 0.058
    arrow(ax, (blue_x[0], ydim), (blue_x[1], ydim), NAVY,
          lw=0.75, scale=6.2, zorder=3)
    arrow(ax, (blue_x[1], ydim), (blue_x[0], ydim), NAVY,
          lw=0.75, scale=6.2, zorder=3)
    ax.text(np.mean(blue_x), 0.028, r"$a$", ha="center", va="center",
            fontsize=7.7, color=INK)


def channel_map(ax, cancelled=False):
    """Floquet channel topology at the verified off-Gamma Rayleigh point."""
    ox, oy = 0.50, 0.405
    rx, ry = 0.275, 0.325

    # Upper half of the real radiation cone |kx,n| <= k.
    ax.add_patch(Arc((ox, oy), 2 * rx, 2 * ry, theta1=0, theta2=180,
                     color=LIGHT, lw=1.05, zorder=0))
    ax.plot([ox - rx - 0.045, ox + rx + 0.085], [oy, oy], color=INK, lw=0.85, zorder=1)
    ax.plot([ox, ox], [oy, oy + ry + 0.035], color=INK, lw=0.85, zorder=1)
    arrow(ax, (ox + rx + 0.055, oy), (ox + rx + 0.085, oy), INK, lw=0.85, scale=6.5)
    arrow(ax, (ox, oy + ry + 0.005), (ox, oy + ry + 0.035), INK, lw=0.85, scale=6.5)
    ax.text(ox + rx + 0.096, oy - 0.003, r"$k_x$", ha="left", va="center", fontsize=7.5)
    ax.text(ox - 0.025, oy + ry + 0.006, r"$k_y$", ha="right", va="center", fontsize=7.6)

    # Normalized channel coordinates: kx,n/k and ky,n/k.
    u0 = KAPPA / OMEGA
    v0 = np.sqrt(1.0 - u0 ** 2)
    p0 = (ox + rx * u0, oy + ry * v0)
    pm = (ox - rx, oy)  # n=-1 lies exactly at the left branch point.
    up = (KAPPA + 1.0) / OMEGA
    pp = (ox + rx * up, oy)  # outside the real cone: ky,+1 is imaginary.

    active_alpha = 0.24 if cancelled else 1.0
    wave_arrow(ax, (ox, oy), p0, BLUE, alpha=active_alpha)
    wave_arrow(ax, (ox, oy), pm, ORANGE, alpha=active_alpha)

    ax.add_patch(Circle((ox, oy), 0.010, facecolor=INK, edgecolor="none", zorder=10))
    ax.add_patch(Circle(pp, 0.010, facecolor=MID, edgecolor="white", lw=0.4, zorder=5))

    # The true n=0 direction is theta=7.10 degrees from the surface normal.
    ax.add_patch(Arc((ox, oy), 0.115, 0.135,
                     theta1=90 - THETA_DEG, theta2=90,
                     color=MID, lw=0.75, zorder=8))
    ax.text(ox + 0.044, oy + 0.105, r"$7.10^\circ$", ha="left", va="center",
            fontsize=7.3, color=MID)

    ax.text(p0[0] + 0.022, p0[1] - 0.002, r"$n=0$", ha="left", va="center",
            fontsize=7.7, color=BLUE)
    ax.text(pm[0] - 0.010, pm[1] + 0.055, r"$n=-1$", ha="center", va="bottom",
            fontsize=7.7, color=ORANGE)
    ax.text(pm[0] - 0.010, pm[1] + 0.017, r"$k_y=0$", ha="center", va="bottom",
            fontsize=7.3, color=MID)
    ax.text(pp[0] - 0.004, pp[1] + 0.042, r"$n=+1$", ha="center", va="bottom",
            fontsize=7.2, color=MID)
    ax.text(pp[0] - 0.004, pp[1] - 0.032, "evanescent", ha="center", va="top",
            fontsize=6.4, color=MID)

    if cancelled:
        zero_mark(ax, p0)
        zero_mark(ax, pm)

    periodic_dimer(ax)


def panel_a(ax):
    setup(ax)
    header(ax, "a", "Conventional BIC")
    localized_mode(ax, center=(0.50, 0.48))

    # Two independent radiation contributions close at one continuum port.
    arrow(ax, (0.415, 0.52), (0.49, 0.745), BLUE, lw=1.45)
    arrow(ax, (0.585, 0.52), (0.51, 0.745), ORANGE, lw=1.45)
    zero_mark(ax, (0.50, 0.755))
    for w, h in [(0.68, 0.53), (0.82, 0.65)]:
        ax.add_patch(Arc((0.50, 0.48), w, h, theta1=22, theta2=158,
                         color=LIGHT, lw=0.75, zorder=0))
        ax.add_patch(Arc((0.50, 0.48), w, h, theta1=202, theta2=338,
                         color=LIGHT, lw=0.75, zorder=0))
    ax.text(0.50, 0.105, r"$A_{\rm rad}=0$", ha="center", va="center",
            fontsize=8.8, color=INK)


def panel_b(ax):
    setup(ax)
    header(ax, "b", "Rayleigh anomaly", r"off-$\Gamma$")
    channel_map(ax, cancelled=False)
    ax.text(0.50, 0.835, r"$A_0\ne0,\quad A_{-1}\ne0$", ha="center", va="center",
            fontsize=8.6, color=INK)


def panel_c(ax):
    setup(ax)
    header(ax, "c", "Rayleigh BIC", r"off-$\Gamma$")
    channel_map(ax, cancelled=True)
    ax.text(0.50, 0.835, r"$A_0=0,\quad A_{-1}=0$", ha="center", va="center",
            fontsize=8.6, color=INK)


def main():
    fig = plt.figure(figsize=(7.2047, 2.92), facecolor="white")
    gs = fig.add_gridspec(1, 3, width_ratios=[0.90, 1.05, 1.05], wspace=0.055)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    panel_a(axes[0])
    panel_b(axes[1])
    panel_c(axes[2])
    fig.subplots_adjust(left=0.012, right=0.992, bottom=0.02, top=0.99)

    out = Path(__file__).resolve().parents[1] / "figures"
    out.mkdir(parents=True, exist_ok=True)
    stem = out / "fig1_potential_corrected_nature"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.015)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.015)
    fig.savefig(stem.with_suffix(".png"), dpi=450, bbox_inches="tight", pad_inches=0.015)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pad_inches=0.015)
    plt.close(fig)


if __name__ == "__main__":
    main()
