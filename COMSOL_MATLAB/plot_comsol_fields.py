#!/usr/bin/env python3
"""Plot the solved COMSOL eigenmode and driven field at final journal size."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "COMSOL_MATLAB" / "rounded_180k_results"
OUT = DATA / "figure_comsol_fields_rounded_180k"
FIGURE_WIDTH_MM = 183
FIGURE_HEIGHT_MM = 78
RASTER_DPI = 600

# Figure contract: the rounded, manufacturable unit cell retains a localized
# high-Q eigenmode and produces the corresponding resonant driven field near
# 180 kHz.  The two panels are complementary eigenproblem and scattering
# evidence, not independent claims.
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.0,
        "axes.labelsize": 7.0,
        "xtick.labelsize": 6.2,
        "ytick.labelsize": 6.2,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)

PRESSURE_CMAP = LinearSegmentedColormap.from_list(
    "pressure_diverging", ["#174A7E", "#6EA6C9", "#F7F7F5", "#F0A05A", "#B6422E"]
)


def load_field(stem: str):
    x = np.loadtxt(DATA / f"{stem}_x_m.csv", delimiter=",") * 1e3
    y = np.loadtxt(DATA / f"{stem}_y_m.csv", delimiter=",") * 1e3
    real = np.loadtxt(DATA / f"{stem}_real.csv", delimiter=",")
    imag = np.loadtxt(DATA / f"{stem}_imag.csv", delimiter=",")
    assert real.shape == (y.size, x.size)
    assert imag.shape == real.shape
    return x, y, real + 1j * imag


def phase_aligned_real(field: np.ndarray) -> np.ndarray:
    """Choose one global phase so the largest finite sample is positive real."""
    finite = np.isfinite(field)
    idx = np.nanargmax(np.abs(field))
    phase = np.angle(field.flat[idx])
    return np.real(field * np.exp(-1j * phase))


def repeat_bloch(x, field, period_mm, kx_rad_m, copies=5):
    shifts = np.arange(-(copies // 2), copies // 2 + 1)
    x_blocks, f_blocks = [], []
    for block_index, shift in enumerate(shifts):
        phase = np.exp(1j * kx_rad_m * (shift * period_mm * 1e-3))
        # Adjacent cells share one boundary coordinate.  Drop the repeated
        # column so the plotting grid remains strictly increasing.
        start = 0 if block_index == 0 else 1
        x_blocks.append(x[start:] + shift * period_mm)
        f_blocks.append(field[:, start:] * phase)
    return np.concatenate(x_blocks), np.concatenate(f_blocks, axis=1)


def add_steel(ax, period_mm, copies=5):
    w1, d1, w2, d2 = 4.41, 5.15, 1.00, 1.32
    x1, x2 = 0.505, 5.915
    for shift in np.arange(-(copies // 2), copies // 2 + 1):
        x0 = shift * period_mm
        ax.add_patch(
            Rectangle((x0, -10.0), period_mm, 10.0, facecolor="#B6BDC3", edgecolor="none", zorder=-3)
        )
        # White openings are subsequently filled by the COMSOL pressure map.
        ax.add_patch(Rectangle((x0 + x1, -d1), w1, d1, facecolor="white", edgecolor="none", zorder=-2))
        ax.add_patch(Rectangle((x0 + x2, -d2), w2, d2, facecolor="white", edgecolor="none", zorder=-2))
        ax.plot([x0, x0 + x1], [0, 0], color="#525B63", lw=0.65, zorder=4)
        ax.plot([x0 + x1 + w1, x0 + x2], [0, 0], color="#525B63", lw=0.65, zorder=4)
        ax.plot([x0 + x2 + w2, x0 + period_mm], [0, 0], color="#525B63", lw=0.65, zorder=4)


def draw_panel(ax, x, y, field, period_mm, kx_rad_m, title, subtitle):
    xr, fr = repeat_bloch(x, field, period_mm, kx_rad_m)
    plotted = phase_aligned_real(fr)
    finite = plotted[np.isfinite(plotted)]
    # Sharp ideal groove corners create isolated FEM pressure spikes.  A
    # symmetric 99.5-percentile display limit preserves the raw CSV while
    # preventing a few corner pixels from erasing the spatial field pattern.
    limit = np.percentile(np.abs(finite), 99.5)
    shown = np.clip(plotted / limit, -1.0, 1.0)
    add_steel(ax, period_mm)
    mesh = ax.pcolormesh(xr, y, shown, cmap=PRESSURE_CMAP, vmin=-1, vmax=1, shading="auto", rasterized=True)
    ax.set_xlim(-2.5 * period_mm, 2.5 * period_mm)
    ax.set_ylim(-6.2, 16.0)
    ax.set_aspect("equal")
    ax.set_title(title, loc="left", fontsize=7.4, fontweight="bold", pad=3.0)
    ax.text(0.0, 14.45, subtitle, ha="center", va="top", fontsize=6.2, color="#33383D")
    ax.set_xlabel("x (mm)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#4A5157")
    ax.tick_params(direction="out", pad=1.5)
    return mesh


def main():
    x_e, y_e, eigen = load_field("eigen_field")
    x_s, y_s, scattering = load_field("scattering_field")
    assert np.all(np.diff(x_e) > 0) and np.all(np.diff(y_e) > 0)
    assert np.all(np.diff(x_s) > 0) and np.all(np.diff(y_s) > 0)
    assert np.allclose(x_e, x_s) and np.allclose(y_e, y_s)

    period_mm = 7.42
    kx_rad_m = 93.1933179121876
    fig, axes = plt.subplots(
        1, 2,
        figsize=(FIGURE_WIDTH_MM / 25.4, FIGURE_HEIGHT_MM / 25.4),
        constrained_layout=False,
    )
    fig.subplots_adjust(left=0.062, right=0.972, bottom=0.18, top=0.89, wspace=0.14)

    mappable = draw_panel(
        axes[0], x_e, y_e, eigen, period_mm, kx_rad_m,
        "a   High-Q eigenmode", "180.046 kHz   |   Q = 252,000",
    )
    draw_panel(
        axes[1], x_s, y_s, scattering, period_mm, kx_rad_m,
        "b   Driven response", "180.040 kHz   |   θ = 7.10°",
    )
    axes[0].set_ylabel("y (mm)")
    axes[1].set_yticklabels([])
    axes[1].set_ylabel("")

    cax = fig.add_axes([0.405, 0.075, 0.19, 0.022])
    cb = fig.colorbar(mappable, cax=cax, orientation="horizontal", ticks=[-1, 0, 1])
    cb.set_label("Phase-aligned Re(p), normalized", labelpad=1.0)
    cb.outline.set_linewidth(0.55)
    cb.ax.tick_params(length=2.0, width=0.5, pad=1.2)

    for suffix, kwargs in {
        ".pdf": {},
        ".svg": {},
        ".png": {"dpi": 300},
        ".tiff": {"dpi": RASTER_DPI},
    }.items():
        fig.savefig(OUT.with_suffix(suffix), bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    main()
