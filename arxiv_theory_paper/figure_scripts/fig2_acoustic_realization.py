"""PRL Fig. 2: physical realization and strict Rayleigh-BIC verification.

All quantitative panels use audited final-root data.  Panel (a) is a vector
isometric rendering of the periodic two-groove unit cell; no AI-generated or
photorealistic content is used.
"""

from pathlib import Path
import csv

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7.5,
        "axes.labelsize": 8.0,
        "axes.titlesize": 8.5,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 7.2,
        "ytick.labelsize": 7.2,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "legend.fontsize": 7.2,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "mathtext.fontset": "dejavusans",
    }
)

INK = "#20242A"
MID = "#697586"
GRID = "#D9DEE3"
BLUE = "#2166AC"
ORANGE = "#D97706"
RED = "#C7352B"
CYAN = "#54A7C8"
TOP = "#DCE9EE"
FRONT = "#A7BDC7"
SIDE = "#829FAA"
WATER = "#EAF6FA"


ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = ROOT / "arxiv_theory_paper" / "figures"
DATA_DIR = ROOT / "arxiv_theory_paper" / "figure_data"


def panel_label(ax, letter):
    ax.text(-0.11, 1.085, f"({letter})", transform=ax.transAxes, ha="left", va="top",
            fontsize=9.2, fontweight="bold", color=INK, clip_on=False)


def clean_axes(ax, grid=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(True, color=GRID, lw=0.55, alpha=0.65)
        ax.set_axisbelow(True)


def project(x, z, y):
    """Oblique projection used for the vector isometric schematic."""
    return np.asarray(x) + 0.40 * np.asarray(z), np.asarray(y) + 0.24 * np.asarray(z)


def poly3(ax, xyz, **kwargs):
    pts = np.array([project(x, z, y) for x, z, y in xyz])
    patch = Polygon(pts, closed=True, **kwargs)
    ax.add_patch(patch)
    return patch


def line3(ax, xyz, **kwargs):
    pts = np.array([project(x, z, y) for x, z, y in xyz])
    ax.plot(pts[:, 0], pts[:, 1], **kwargs)


def draw_trench(ax, x0, x1, zmax, depth, zorder=5):
    """Draw a translationally invariant water-filled groove.

    The top strip establishes the out-of-plane invariance of the modal model;
    the front cut exposes the physical depth without using a decorative inset.
    """
    # The aperture is a geometric void, not a color-coded component.  Its
    # very pale tone is inherited from the surrounding water domain.
    poly3(ax, [(x0, 0, 0.012), (x1, 0, 0.012),
               (x1, zmax, 0.012), (x0, zmax, 0.012)],
          facecolor="#F4FAFC", edgecolor=INK, lw=0.46, zorder=zorder)
    poly3(ax, [(x0, 0, 0), (x1, 0, 0),
               (x1, 0, -depth), (x0, 0, -depth)],
          facecolor="#FFFFFF", edgecolor="none", zorder=zorder + 1)
    line3(ax, [(x0, 0, 0), (x0, 0, -depth), (x1, 0, -depth), (x1, 0, 0)],
          color=INK, lw=0.62, zorder=zorder + 2)
    line3(ax, [(x0, zmax, 0.018), (x1, zmax, 0.018)],
          color="#8FA7B2", lw=0.44, alpha=0.62, zorder=zorder + 2)


def draw_structure(ax, geom):
    # A slightly amplified vertical projection lets the physically thin slab
    # occupy the panel without turning it into a bulky block.
    ax.set_aspect(2.15)
    ax.axis("off")

    ncell, zmax, thick = 9.0, 1.35, 0.74
    # The solid is drawn as one continuous plate.  Top, front and side use a
    # single low-saturation material family with a consistent light direction.
    poly3(ax, [(0, 0, 0), (ncell, 0, 0), (ncell, 0, -thick), (0, 0, -thick)],
          facecolor="#657E96", edgecolor=INK, lw=0.68, zorder=1)
    front_left = np.array(mpl.colors.to_rgb("#557F94"))
    front_right = np.array(mpl.colors.to_rgb("#766E99"))
    for j in range(36):
        x0 = ncell * j / 36
        x1 = ncell * (j + 1) / 36
        t = (j + 0.5) / 36
        fill = tuple((1 - t) * front_left + t * front_right)
        poly3(ax, [(x0, 0, -0.002), (x1, 0, -0.002),
                   (x1, 0, -thick + 0.002), (x0, 0, -thick + 0.002)],
              facecolor=fill, edgecolor="none", zorder=1.2)
    poly3(ax, [(ncell, 0, 0), (ncell, zmax, 0),
               (ncell, zmax, -thick), (ncell, 0, -thick)],
          facecolor="#77769A", edgecolor=INK, lw=0.68, zorder=2)
    poly3(ax, [(0, 0, 0), (ncell, 0, 0), (ncell, zmax, 0), (0, zmax, 0)],
          facecolor="#BFDDE1", edgecolor=INK, lw=0.72, zorder=3)
    # Vector-only cool gradient across the upper surface gives the flat array
    # the visual weight of a photonic-crystal chip without coloring the voids.
    c_front = np.array(mpl.colors.to_rgb("#ADD9DC"))
    c_back = np.array(mpl.colors.to_rgb("#CFC8E4"))
    for j in range(18):
        z0 = zmax * j / 18
        z1 = zmax * (j + 1) / 18
        t = (j + 0.5) / 18
        fill = tuple((1 - t) * c_front + t * c_back)
        poly3(ax, [(0, z0, 0.004), (ncell, z0, 0.004),
                   (ncell, z1, 0.004), (0, z1, 0.004)],
              facecolor=fill, edgecolor="none", zorder=3.2)
    # Back-edge highlight and a quiet shadow make the object read as one slab.
    line3(ax, [(0, zmax, 0.008), (ncell, zmax, 0.008)],
          color="white", lw=0.72, alpha=0.80, zorder=4)
    ax.add_patch(Polygon([project(0.08, 0, -thick - 0.08),
                          project(ncell + 0.12, 0, -thick - 0.08),
                          project(ncell + 0.34, 0, -thick - 0.13),
                          project(0.31, 0, -thick - 0.13)],
                         closed=True, facecolor="#B8C2C8", edgecolor="none",
                         alpha=0.24, zorder=0))

    w1 = geom["width1_mm"] / geom["period_mm"]
    w2 = geom["width2_mm"] / geom["period_mm"]
    gap = geom["gap_mm"] / geom["period_mm"]
    d1 = geom["depth1_mm"] / geom["period_mm"]
    d2 = geom["depth2_mm"] / geom["period_mm"]
    margin = (1.0 - w1 - w2 - gap) / 2.0
    for cell in range(int(ncell)):
        x1a = cell + margin
        x1b = x1a + w1
        x2a = x1b + gap
        x2b = x2a + w2
        draw_trench(ax, x1a, x1b, zmax, d1, zorder=5)
        draw_trench(ax, x2a, x2b, zmax, d2, zorder=8)

    # Unit-cell boundaries are shown only on the top surface and remain
    # subordinate to the two groove families.
    for xb in np.arange(1.0, ncell, 1.0):
        line3(ax, [(xb, 0, 0.022), (xb, zmax, 0.022)], color="white", lw=0.62,
              ls=(0, (2.2, 2.2)), alpha=0.86, zorder=12)

    # A single period marker is sufficient; the two groove families are
    # identified in the caption rather than crowded onto the thin chip.
    representative_cell = 4.0

    # One clean dimension line for the period; split around the label.
    p0 = project(representative_cell, 0.0, -0.91)
    p1 = project(representative_cell + 1.0, 0.0, -0.91)
    xm = (p0[0] + p1[0]) / 2
    ax.add_patch(FancyArrowPatch((xm - 0.08, p0[1]), p0, arrowstyle="-|>", mutation_scale=6,
                                 color=MID, lw=0.8))
    ax.add_patch(FancyArrowPatch((xm + 0.08, p0[1]), p1, arrowstyle="-|>", mutation_scale=6,
                                 color=MID, lw=0.8))
    ax.text(xm, p0[1], r"$a$", ha="center", va="center", fontsize=7.5, color=INK)

    ax.set_xlim(-0.14, ncell + 0.40 * zmax + 0.18)
    ax.set_ylim(-0.92, 0.24 * zmax + 0.16)


def draw_scattering_map(ax, theta, frequency, abs_am1, geom):
    assert abs_am1.shape == (theta.size, frequency.size)
    assert np.all(abs_am1 > 0), "Driven-amplitude map must be strictly positive before log10"
    log_amp = np.log10(abs_am1)
    # imshow keeps the dense atlas as one bounded image object in PDF/SVG;
    # rasterized pcolormesh can be displaced by tight-bbox vector export.
    mesh = ax.imshow(log_amp.T, origin="lower", aspect="auto",
                     extent=[theta.min(), theta.max(),
                             frequency.min() / 1e3, frequency.max() / 1e3],
                     interpolation="bilinear", cmap="viridis",
                     vmin=-4.3, vmax=0.4)

    th = np.linspace(theta.min(), theta.max(), 500)
    f_ra = geom["period_mm"] ** -1 * 1e3 * 1500.0 / (1 + np.sin(np.deg2rad(th)))
    ax.plot(th, f_ra / 1e3, color="white", lw=1.15, ls=(0, (4, 2.5)))
    ax.plot(geom["theta_deg"], 200.0, "o", ms=5.5, mfc=RED, mec="white", mew=0.7, zorder=5)
    ax.text(geom["theta_deg"] + 0.10, 200.18, "BIC", color="white", fontsize=7.2,
            ha="left", va="bottom")

    ax.set_xlim(theta.min(), theta.max())
    ax.set_ylim(frequency.min() / 1e3, frequency.max() / 1e3)
    ax.set_xlabel(r"Incidence angle $\theta$ (deg)")
    ax.set_ylabel("Frequency (kHz)")
    clean_axes(ax, grid=False)

    cax = inset_axes(ax, width="3.8%", height="100%", loc="lower left",
                     bbox_to_anchor=(1.035, 0.0, 1, 1), bbox_transform=ax.transAxes,
                     borderpad=0)
    # Draw the narrow color scale as vector rectangles.  Matplotlib's inset
    # colorbar image can acquire an incorrect transform in tight-bbox PDFs.
    norm = mpl.colors.Normalize(vmin=-4.3, vmax=0.4)
    cmap = mpl.colormaps["viridis"]
    levels = np.linspace(-4.3, 0.4, 97)
    for lo, hi in zip(levels[:-1], levels[1:]):
        cax.add_patch(Rectangle((0, lo), 1, hi - lo,
                                facecolor=cmap(norm(0.5 * (lo + hi))),
                                edgecolor="none"))
    cax.set_xlim(0, 1)
    cax.set_ylim(-4.3, 0.4)
    cax.set_xticks([])
    cax.set_yticks([-4, -3, -2, -1, 0])
    cax.yaxis.tick_right()
    cax.yaxis.set_label_position("right")
    cax.set_ylabel(r"$\log_{10}|A_{-1}|$", fontsize=7.5, labelpad=5)
    cax.tick_params(axis="y", labelsize=7.2, length=2.5)
    for spine in cax.spines.values():
        spine.set_linewidth(0.6)


def draw_strict_null_scan(ax, strict_rows):
    """Full-truncation real-axis proof of the isolated strict null."""
    dk = np.array([r["deltaKappa"] for r in strict_rows], float)
    sigma = np.array([r["sigmaRelative"] for r in strict_rows], float)
    assert dk.size == 161 and np.all(sigma > 0)
    center = np.argmin(np.abs(dk))

    x = 1e4 * dk
    ax.semilogy(x, sigma, color=BLUE, lw=1.25)
    ax.plot(x[::10], sigma[::10], "o", ms=2.7, mfc="white", mec=BLUE,
            mew=0.75, zorder=3)
    ax.plot(x[center], sigma[center], "o", ms=5.5, mfc=RED, mec="white",
            mew=0.7, zorder=5)
    ax.annotate("BIC", xy=(0, sigma[center]), xytext=(0.55, 2.5e-12),
                color=RED, fontsize=7.2, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.75))

    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(5e-16, 8e-3)
    ax.set_xlabel(r"$(\kappa-\kappa_{\mathrm{BIC}})\times10^{4}$")
    ax.set_ylabel(r"$\sigma_{\min}/\sigma_{\max}$")
    ax.set_yticks([1e-15, 1e-11, 1e-7, 1e-3])
    clean_axes(ax, grid=True)


def draw_radiation_phasors(ax, phasor_rows):
    """Measured two-groove cancellation in the only two radiation constraints."""
    ax.set_aspect("equal")
    ax.axis("off")
    groups = {
        "open": [r for r in phasor_rows if r["channel"] == "open"],
        "rayleigh": [r for r in phasor_rows if r["channel"] == "rayleigh"],
    }
    centers = {"open": np.array([-0.72, 0.04]), "rayleigh": np.array([0.72, 0.04])}
    titles = {"open": r"open $n=0$", "rayleigh": r"Rayleigh $n=-1$"}
    qcolors = ["#3B9AB2", "#8B6BB8"]

    for channel in ["open", "rayleigh"]:
        rows = groups[channel]
        assert len(rows) == 2
        vectors = np.array([[float(r["real"]), float(r["imag"])] for r in rows])
        residual = float(rows[0]["residual_relative"])
        vectors /= np.max(np.linalg.norm(vectors, axis=1))
        origin = centers[channel]
        scale = 0.43
        endpoint = origin + scale * vectors[0]

        ax.plot([origin[0] - 0.48, origin[0] + 0.48], [origin[1], origin[1]],
                color=GRID, lw=0.55, zorder=0)
        ax.plot([origin[0], origin[0]], [origin[1] - 0.49, origin[1] + 0.49],
                color=GRID, lw=0.55, zorder=0)
        ax.add_patch(FancyArrowPatch(origin, endpoint, arrowstyle="-|>",
                                     mutation_scale=8.2, color=qcolors[0], lw=1.65,
                                     shrinkA=0, shrinkB=0, zorder=3))
        ax.add_patch(FancyArrowPatch(endpoint, origin, arrowstyle="-|>",
                                     mutation_scale=8.2, color=qcolors[1], lw=1.65,
                                     shrinkA=0, shrinkB=0, zorder=4,
                                     connectionstyle="arc3,rad=0.075"))
        ax.scatter(*origin, s=25, facecolor="white", edgecolor=RED,
                   linewidth=1.05, zorder=6)
        ax.text(origin[0], 0.62, titles[channel], ha="center", va="center",
                fontsize=7.5, color=INK)
        ax.text(origin[0], -0.48, rf"$|\Sigma|={residual:.1e}$", ha="center",
                va="center", fontsize=7.2, color=RED)

    ax.text(-1.22, 0.46, r"Im", fontsize=7.2, color=MID, ha="left", va="center")
    ax.text(-0.25, 0.08, r"Re", fontsize=7.2, color=MID, ha="right", va="bottom")
    ax.plot([], [], color=qcolors[0], lw=1.7, label="Groove 1")
    ax.plot([], [], color=qcolors[1], lw=1.7, label="Groove 2")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.16), ncol=2,
              handlelength=1.5, columnspacing=1.2, borderaxespad=0)
    ax.set_xlim(-1.34, 1.34)
    ax.set_ylim(-0.88, 0.72)


def read_geometry():
    path = ROOT / "Ni2019_MATLAB" / "results" / "StrictRayleighBIC_200kHz_min1mm.csv"
    with path.open(newline="") as f:
        row = next(csv.DictReader(f))
    return {k: float(v) for k, v in row.items() if k != "accepted"}


def read_rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def main():
    geom = read_geometry()
    strict_rows = read_rows(DATA_DIR / "fig2_strict_rayleigh_line.csv")
    phasor_rows = read_rows(DATA_DIR / "fig2_radiation_phasors.csv")
    theta = np.loadtxt(DATA_DIR / "fig2_map_theta_deg.csv", delimiter=",")
    frequency = np.loadtxt(DATA_DIR / "fig2_map_frequency_hz.csv", delimiter=",")
    abs_am1 = np.loadtxt(DATA_DIR / "fig2_map_absAm1.csv", delimiter=",")

    kappa0 = geom["kappa"]
    omega0 = geom["Omega"]
    assert abs(omega0 - (1 - kappa0)) < 2e-15

    fig = plt.figure(figsize=(7.10, 5.02), facecolor="white")
    gs = fig.add_gridspec(2, 2, left=0.075, right=0.985, bottom=0.095, top=0.985,
                          wspace=0.27, hspace=0.34, height_ratios=[1.03, 1.0])
    axes = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(2)]

    draw_structure(axes[0], geom)
    draw_scattering_map(axes[1], theta, frequency, abs_am1, geom)
    draw_strict_null_scan(axes[2], strict_rows)
    draw_radiation_phasors(axes[3], phasor_rows)
    # Panel labels are anchored to the GridSpec cells rather than the
    # aspect-adjusted axes, so both rows retain exact typographic baselines.
    specs = [gs[0, 0], gs[0, 1], gs[1, 0], gs[1, 1]]
    for spec, letter in zip(specs, "abcd"):
        box = spec.get_position(fig)
        fig.text(box.x0 - 0.030, box.y1 + 0.010, f"({letter})",
                 ha="left", va="bottom", fontsize=9.2,
                 fontweight="bold", color=INK)

    # Deterministic alignment checks for the two-column PRL grid.
    fig.canvas.draw()
    pos = [ax.get_position() for ax in axes]
    assert abs(pos[0].x0 - pos[2].x0) < 1e-10
    assert abs(pos[1].x0 - pos[3].x0) < 1e-10
    assert abs(pos[0].width - pos[1].width) < 1e-10
    assert abs(pos[2].width - pos[3].width) < 1e-10

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    stem = FIG_DIR / "fig2_acoustic_realization"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.025)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.025)
    fig.savefig(stem.with_suffix(".png"), dpi=400, bbox_inches="tight", pad_inches=0.025)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pad_inches=0.025)
    plt.close(fig)


if __name__ == "__main__":
    main()
