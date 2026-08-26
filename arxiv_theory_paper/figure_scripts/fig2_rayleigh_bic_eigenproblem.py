#!/usr/bin/env python3
"""Fig. 2 | Realization and eigenmode formation of the acoustic Rayleigh BIC.

The numerical source is MATLAB modal matching.  This Python file performs
publication rendering only.  Its evidence chain is:

  2D z-invariant geometry -> sheet-resolved pole bands -> Q divergence at
  the n=-1 Rayleigh crossing -> simultaneous A0 and A-1 cancellation ->
  homogeneous eigenfield evolution.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "Ni2019_MATLAB" / "results"
FIG_DIR = ROOT / "arxiv_theory_paper" / "figures"
GEOMETRY_MAT = DATA / "StrictRayleighBIC_180kHz_7p10deg_final.mat"
EVOLUTION_CSV = (
    DATA
    / "kx_channel_evolution_physical_180k"
    / "physical_sheet_channel_evolution.csv"
)
TARGET_BZ_CSV = (
    DATA
    / "target_branch_physical_bz_180k"
    / "target_branch_physical_bz.csv"
)
PHYSICAL_POLE_POINTS_CSV = (
    DATA
    / "full_complex_band_180k"
    / "physical_outgoing_pole_points.csv"
)
PHYSICAL_POLE_SEGMENTS_CSV = (
    DATA
    / "full_complex_band_180k"
    / "physical_outgoing_pole_segments.csv"
)
GUIDED_NEAR_GAMMA_CSV = (
    DATA
    / "full_complex_band_180k"
    / "low_guided_branch_near_gamma.csv"
)
FIELD_MAT = (
    DATA
    / "fig2_three_state_fields_180k"
    / "eigenmode_fields_gamma_bic_beyond.mat"
)

MM = 1.0 / 25.4
WIDTH_MM = 178.0
HEIGHT_MM = 205.0

INK = "#20262C"
MUTED = "#68747C"
GRID = "#DDE2E5"
SLAB = "#AFC3C8"
SLAB_TOP = "#C4D4D7"
SLAB_FRONT = "#718A94"
SLAB_SIDE = "#879FA7"
NAVY = "#174E7A"
BLUE = "#2878A8"
TEAL = "#147D78"
PLUM = "#8A4F7D"
ORANGE = "#D77A1F"
RED = "#D43F34"
WHITE = "#FFFFFF"

Q_MAX = 12.0
Q_CMAP = LinearSegmentedColormap.from_list(
    "quiet_to_radiant_q",
    ["#C9D0D4", "#80AFC1", "#278F98", "#F0B23A", "#DC4638"],
)
Q_NORM = Normalize(vmin=0.0, vmax=Q_MAX, clip=True)

FIELD_CMAP = LinearSegmentedColormap.from_list(
    "field_navy_teal_gold",
    ["#111D35", "#174A6A", "#16827D", "#72B776", "#F0CE64"],
)
FIELD_CMAP.set_bad("#B8C5CA")
CHANNEL_REAL_CMAP = LinearSegmentedColormap.from_list(
    "channel_blue_white_orange",
    ["#173F73", "#6EA6C9", "#F7F7F5", "#E9A05B", "#A93632"],
)

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7.2,
        "axes.labelsize": 7.5,
        "axes.titlesize": 7.5,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.2,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 2.7,
        "ytick.major.size": 2.7,
        "xtick.major.width": 0.65,
        "ytick.major.width": 0.65,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "mathtext.fontset": "dejavusans",
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
)


def panel_label(ax: plt.Axes, letter: str, x: float = -0.105, y: float = 1.05) -> None:
    ax.text(
        x,
        y,
        f"({letter})",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.1,
        fontweight="bold",
        color=INK,
        clip_on=False,
        zorder=50,
    )


def style_axis(ax: plt.Axes, grid: bool = False) -> None:
    ax.spines["left"].set_color(INK)
    ax.spines["bottom"].set_color(INK)
    ax.tick_params(colors=INK, pad=2.0)
    if grid:
        ax.grid(color=GRID, linewidth=0.45, alpha=0.75)
        ax.set_axisbelow(True)


def add_continuous_q_curve(
    ax: plt.Axes, branch: list[dict[str, float]], linewidth: float
) -> None:
    """Draw one verified continuous branch despite nonuniform k sampling."""
    ordered = sorted(branch, key=lambda row: row["kappa"])
    x = np.asarray([row["kappa"] for row in ordered], dtype=float)
    y = np.asarray([row["Omega_real"] for row in ordered], dtype=float)
    q = np.asarray([row["log10_Q_display"] for row in ordered], dtype=float)
    assert len(x) >= 2 and np.all(np.diff(x) > 0)
    points = np.column_stack([x, y])
    segments = np.stack([points[:-1], points[1:]], axis=1)
    collection = LineCollection(
        segments,
        cmap=Q_CMAP,
        norm=Q_NORM,
        linewidth=linewidth,
        capstyle="round",
        joinstyle="round",
        zorder=4.2,
    )
    collection.set_array(0.5 * (q[:-1] + q[1:]))
    ax.add_collection(collection)


def load_physical_pole_spectrum() -> tuple[np.ndarray, np.ndarray]:
    """Load independently solved physical-boundary poles and safe local joins."""
    points = np.genfromtxt(PHYSICAL_POLE_POINTS_CSV, delimiter=",", names=True)
    segments = np.genfromtxt(PHYSICAL_POLE_SEGMENTS_CSV, delimiter=",", names=True)
    assert np.nanmin(segments["mode_overlap"]) >= 0.95
    assert np.nanmax(segments["delta_Omega"]) <= 0.035
    return points, segments


def _is_target_fragment(kappa: np.ndarray, omega: np.ndarray) -> np.ndarray:
    """Identify the low-resolution census copy of the dedicated target branch."""
    return (kappa >= 0.30) & (omega >= 0.88) & (omega <= 0.91)


def add_physical_pole_segments(
    ax: plt.Axes, points: np.ndarray, segments: np.ndarray
) -> None:
    """Render only overlap- and sheet-consistent local pole connections."""
    target_duplicate = _is_target_fragment(
        segments["kappa_0"], segments["Omega_real_0"]
    ) & _is_target_fragment(segments["kappa_1"], segments["Omega_real_1"])
    keep = ~target_duplicate
    s = segments[keep]
    radiative = s["sheet_signature"] > 0
    guided = ~radiative

    for sign in (-1.0, 1.0):
        if np.any(guided):
            x0 = sign * s["kappa_0"][guided]
            x1 = sign * s["kappa_1"][guided]
            grey_segments = np.stack(
                [
                    np.column_stack([x0, s["Omega_real_0"][guided]]),
                    np.column_stack([x1, s["Omega_real_1"][guided]]),
                ],
                axis=1,
            )
            ax.add_collection(
                LineCollection(
                    grey_segments,
                    colors="#727E86",
                    linewidths=0.85,
                    alpha=0.82,
                    capstyle="round",
                    zorder=2.8,
                )
            )
        if np.any(radiative):
            x0 = sign * s["kappa_0"][radiative]
            x1 = sign * s["kappa_1"][radiative]
            coloured_segments = np.stack(
                [
                    np.column_stack([x0, s["Omega_real_0"][radiative]]),
                    np.column_stack([x1, s["Omega_real_1"][radiative]]),
                ],
                axis=1,
            )
            q = 0.5 * (
                np.log10(np.clip(s["Q_0"][radiative], 1.0, np.inf))
                + np.log10(np.clip(s["Q_1"][radiative], 1.0, np.inf))
            )
            lc = LineCollection(
                coloured_segments,
                cmap=Q_CMAP,
                norm=Q_NORM,
                linewidths=0.82,
                alpha=0.78,
                capstyle="round",
                zorder=3.0,
            )
            lc.set_array(np.clip(q, 0.0, Q_MAX))
            ax.add_collection(lc)

    point_keep = ~_is_target_fragment(points["kappa"], points["Omega_real"])
    p = points[point_keep]
    for sign in (-1.0, 1.0):
        nonzero = np.ones(len(p), dtype=bool) if sign > 0 else p["kappa"] > 1e-12
        pg = p[nonzero & (p["guided"] > 0.5)]
        pr = p[nonzero & (p["guided"] < 0.5)]
        if len(pg):
            ax.scatter(
                sign * pg["kappa"],
                pg["Omega_real"],
                s=4.6,
                color="#727E86",
                linewidths=0,
                alpha=0.88,
                zorder=3.2,
            )
        if len(pr):
            q = np.clip(
                np.log10(np.clip(pr["Q"], 1.0, np.inf)), 0.0, Q_MAX
            )
            ax.scatter(
                sign * pr["kappa"],
                pr["Omega_real"],
                c=q,
                cmap=Q_CMAP,
                norm=Q_NORM,
                s=5.2,
                linewidths=0,
                alpha=0.92,
                zorder=3.4,
            )

    # A physical-boundary pole ends at a square-root cut.  Open endpoints
    # make this topology explicit and prevent separate segments from being
    # read as a single band that numerically jumps.
    degree: dict[int, int] = {}
    for id0, id1, retained in zip(
        segments["point_id_0"].astype(int),
        segments["point_id_1"].astype(int),
        keep,
    ):
        if not retained:
            continue
        degree[id0] = degree.get(id0, 0) + 1
        degree[id1] = degree.get(id1, 0) + 1
    is_segment_end = np.asarray(
        [degree.get(int(pid), 0) <= 1 for pid in points["point_id"]], dtype=bool
    )
    near_cut = (
        is_segment_end
        & (points["cut_distance"] < 0.012)
        & ~_is_target_fragment(points["kappa"], points["Omega_real"])
        & (points["kappa"] > 0.20)
    )
    endpoint_candidates = points[near_cut]
    if len(endpoint_candidates):
        for sign in (-1.0, 1.0):
            ax.scatter(
                sign * endpoint_candidates["kappa"],
                endpoint_candidates["Omega_real"],
                s=15,
                facecolor=WHITE,
                edgecolor="#68747C",
                linewidth=0.65,
                zorder=8,
            )

    # The lowest guided branch approaches the n=0 continuum at Gamma.  The
    # dedicated continuation resolves its near-threshold asymptote.
    near = np.genfromtxt(GUIDED_NEAR_GAMMA_CSV, delimiter=",", names=True)
    positive_k = np.r_[near["kappa"], 0.04]
    positive_o = np.r_[near["Omega_real"], 0.0398231220045979]
    order = np.argsort(positive_k)
    positive_k = positive_k[order]
    positive_o = positive_o[order]
    for sign in (-1.0, 1.0):
        ax.plot(
            sign * positive_k,
            positive_o,
            color="#727E86",
            linewidth=0.85,
            alpha=0.82,
            solid_capstyle="round",
            zorder=2.8,
        )
    ax.scatter(
        [0.0],
        [0.0],
        s=14,
        facecolor=WHITE,
        edgecolor="#727E86",
        linewidth=0.65,
        zorder=8,
    )


def load_geometry() -> dict[str, float]:
    raw = loadmat(GEOMETRY_MAT, squeeze_me=True, struct_as_record=False)
    cfg = raw["cfg"]
    kappa = float(np.asarray(raw["x"])[0])
    omega = 1.0 - kappa
    a_mm = float(raw["a"]) * 1e3
    widths = np.asarray(cfg.widths, dtype=float)
    depths = np.asarray(cfg.depths, dtype=float)
    gaps = np.atleast_1d(np.asarray(cfg.gaps, dtype=float))
    occupied = float(widths.sum() + gaps.sum())
    margin = 0.5 * (1.0 - occupied)
    geometry = {
        "kappa": kappa,
        "omega": omega,
        "theta_deg": float(np.degrees(np.arcsin(kappa / omega))),
        "a_mm": a_mm,
        "w1": float(widths[0]),
        "w2": float(widths[1]),
        "gap": float(gaps[0]),
        "d1": float(depths[0]),
        "d2": float(depths[1]),
        "margin": margin,
        "w1_mm": float(widths[0] * a_mm),
        "w2_mm": float(widths[1] * a_mm),
        "gap_mm": float(gaps[0] * a_mm),
        "d1_mm": float(depths[0] * a_mm),
        "d2_mm": float(depths[1] * a_mm),
    }
    assert geometry["margin"] > 0
    assert np.isclose(geometry["omega"], 1.0 - geometry["kappa"])
    return geometry


def arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str = INK,
    both: bool = False,
    width: float = 0.65,
    scale: float = 6.5,
    zorder: int = 20,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="<->" if both else "-|>",
            mutation_scale=scale,
            linewidth=width,
            color=color,
            shrinkA=0,
            shrinkB=0,
            zorder=zorder,
        )
    )


def draw_dimension(
    ax: plt.Axes,
    p0: tuple[float, float],
    p1: tuple[float, float],
    label: str,
    text_offset: tuple[float, float] = (0.0, 0.0),
    rotate: float = 0.0,
) -> None:
    arrow(ax, p0, p1, color=MUTED, both=True, width=0.6, scale=5.2)
    xm = 0.5 * (p0[0] + p1[0]) + text_offset[0]
    ym = 0.5 * (p0[1] + p1[1]) + text_offset[1]
    ax.text(
        xm,
        ym,
        label,
        ha="center",
        va="center",
        fontsize=7.2,
        color=INK,
        rotation=rotate,
        rotation_mode="anchor",
        zorder=25,
    )


def draw_geometry(ax: plt.Axes, g: dict[str, float]) -> None:
    """Thin z-invariant extrusion plus a dimensioned x-y unit-cell section."""
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")

    # A deliberately flat extrusion: its small front face is visual depth,
    # not a modeled finite z thickness.  Every aperture is an unfilled void.
    origin = np.array([0.055, 0.682])
    vx = np.array([0.147, 0.0])
    vz = np.array([0.165, 0.128])
    down = np.array([0.0, -0.025])
    n_cells = 5
    top = np.array(
        [origin, origin + n_cells * vx, origin + n_cells * vx + vz, origin + vz]
    )
    front = np.array(
        [origin, origin + n_cells * vx, origin + n_cells * vx + down, origin + down]
    )
    side = np.array(
        [
            origin + n_cells * vx,
            origin + n_cells * vx + vz,
            origin + n_cells * vx + vz + down,
            origin + n_cells * vx + down,
        ]
    )
    ax.add_patch(Polygon(front, facecolor=SLAB_FRONT, edgecolor=INK, linewidth=0.55))
    ax.add_patch(Polygon(side, facecolor=SLAB_SIDE, edgecolor=INK, linewidth=0.55))
    ax.add_patch(Polygon(top, facecolor=SLAB_TOP, edgecolor=INK, linewidth=0.65))

    for cell in range(n_cells):
        cell_origin = origin + cell * vx
        for start, width in (
            (g["margin"], g["w1"]),
            (g["margin"] + g["w1"] + g["gap"], g["w2"]),
        ):
            p = cell_origin + start * vx + 0.075 * vz
            slot = np.array(
                [p, p + width * vx, p + width * vx + 0.84 * vz, p + 0.84 * vz]
            )
            ax.add_patch(
                Polygon(slot, facecolor=WHITE, edgecolor="#51636A", linewidth=0.48)
            )

    # Highlight exactly one period without visually interrupting the array.
    central = 2
    p = origin + central * vx
    cell_outline = np.array([p, p + vx, p + vx + vz, p + vz, p])
    ax.plot(cell_outline[:, 0], cell_outline[:, 1], color=TEAL, linewidth=0.95)
    # Draw continuation dots as vector primitives so the PDF never depends
    # on a font-specific ellipsis glyph.
    ax.scatter(
        [0.010, 0.025, 0.040, 0.960, 0.975, 0.990],
        [0.765] * 6,
        s=3.0,
        color=MUTED,
        linewidths=0,
        clip_on=False,
        zorder=12,
    )

    # The slots continue indefinitely along z.  This is the only 3D cue.
    z0 = origin + 3.55 * vx + 0.70 * vz
    z1 = z0 + 0.42 * vz
    arrow(ax, tuple(z0), tuple(z1), color=INK, width=0.7, scale=7.0)
    ax.text(z1[0] + 0.012, z1[1] + 0.006, r"$z$", fontsize=7.3, color=INK)
    ax.text(
        z1[0] - 0.010,
        z1[1] + 0.048,
        "infinite along z",
        fontsize=7.2,
        color=MUTED,
        ha="center",
    )

    # Period dimension aligned with the highlighted cell.
    yp = 0.625
    x_left = float((origin + central * vx)[0])
    x_right = float((origin + (central + 1) * vx)[0])
    ax.plot([x_left, x_left], [yp + 0.012, origin[1]], color=MUTED, linewidth=0.45)
    ax.plot([x_right, x_right], [yp + 0.012, origin[1]], color=MUTED, linewidth=0.45)
    draw_dimension(
        ax,
        (x_left, yp),
        (x_right, yp),
        f"a = {g['a_mm']:.2f} mm",
        text_offset=(0.0, -0.025),
    )

    # Dimensioned x-y cross-section.  The white rectangles are grooves;
    # the surrounding blue-gray region is the rigid substrate.
    sec_left, sec_right = 0.105, 0.895
    sec_width = sec_right - sec_left
    surface_y = 0.375
    base_y = 0.055
    depth_scale = (surface_y - base_y) / 0.75
    ax.add_patch(
        Rectangle(
            (sec_left, base_y),
            sec_width,
            surface_y - base_y,
            facecolor=SLAB,
            edgecolor=INK,
            linewidth=0.65,
        )
    )
    gx1 = sec_left + sec_width * g["margin"]
    gx2 = sec_left + sec_width * (g["margin"] + g["w1"] + g["gap"])
    gw1, gw2 = sec_width * g["w1"], sec_width * g["w2"]
    gd1, gd2 = depth_scale * g["d1"], depth_scale * g["d2"]
    ax.add_patch(
        Rectangle(
            (gx1, surface_y - gd1),
            gw1,
            gd1,
            facecolor=WHITE,
            edgecolor=INK,
            linewidth=0.6,
        )
    )
    ax.add_patch(
        Rectangle(
            (gx2, surface_y - gd2),
            gw2,
            gd2,
            facecolor=WHITE,
            edgecolor=INK,
            linewidth=0.6,
        )
    )
    ax.plot([sec_left, sec_right], [surface_y, surface_y], color=INK, linewidth=0.72)

    y_width = 0.423
    for xa, xb in ((gx1, gx1 + gw1), (gx2, gx2 + gw2)):
        ax.plot([xa, xa], [surface_y + 0.006, y_width - 0.006], color=MUTED, linewidth=0.4)
        ax.plot([xb, xb], [surface_y + 0.006, y_width - 0.006], color=MUTED, linewidth=0.4)
    draw_dimension(
        ax,
        (gx1, y_width),
        (gx1 + gw1, y_width),
        rf"$w_1={g['w1_mm']:.2f}$",
        text_offset=(0.0, 0.022),
    )
    draw_dimension(
        ax,
        (gx1 + gw1, y_width - 0.041),
        (gx2, y_width - 0.041),
        rf"$g={g['gap_mm']:.2f}$",
        text_offset=(0.0, -0.018),
    )
    draw_dimension(
        ax,
        (gx2, y_width),
        (gx2 + gw2, y_width),
        rf"$w_2={g['w2_mm']:.2f}$",
        text_offset=(0.0, 0.022),
    )

    x_d1 = gx1 - 0.030
    x_d2 = 0.935
    draw_dimension(
        ax,
        (x_d1, surface_y),
        (x_d1, surface_y - gd1),
        rf"$d_1={g['d1_mm']:.2f}$",
        text_offset=(-0.035, 0.0),
        rotate=90,
    )
    draw_dimension(
        ax,
        (x_d2, surface_y),
        (x_d2, surface_y - gd2),
        rf"$d_2={g['d2_mm']:.2f}$",
        text_offset=(0.026, 0.0),
        rotate=90,
    )
    ax.text(
        sec_right,
        base_y - 0.026,
        "dimensions in mm",
        fontsize=6.0,
        color=MUTED,
        ha="right",
        va="top",
    )


def draw_full_band(ax: plt.Axes, cax: plt.Axes, values: np.ndarray) -> None:
    kappa = np.linspace(-0.5, 0.5, 1001)
    ax.fill_between(kappa, 0.0, np.abs(kappa), color="#F1F3F4", linewidth=0)
    ax.plot(
        kappa,
        np.abs(kappa),
        color="#87929B",
        linewidth=0.8,
        linestyle=(0, (3.8, 2.8)),
    )
    ax.plot(
        kappa,
        1.0 - np.abs(kappa),
        color=ORANGE,
        linewidth=1.0,
        linestyle=(0, (4.6, 2.8)),
    )

    pole_points, pole_segments = load_physical_pole_spectrum()
    add_physical_pole_segments(ax, pole_points, pole_segments)

    # Replace the old target/adjacent-sheet continuation by a directly
    # recomputed all-channel outgoing branch. Closed orders use the decaying
    # physical-boundary sign; open orders use the outgoing sign.
    positive_branch = [
        {
            "kappa": float(row["kappa"]),
            "Omega_real": float(row["Omega_real"]),
            "log10_Q_display": float(
                np.clip(
                    np.log10(np.clip(row["Q"], 1.0, np.inf)),
                    0.0,
                    Q_MAX,
                )
            ),
        }
        for row in values
    ]
    negative_branch = [dict(row, kappa=-row["kappa"]) for row in positive_branch]
    add_continuous_q_curve(ax, positive_branch, linewidth=2.15)
    add_continuous_q_curve(ax, negative_branch, linewidth=2.15)

    gamma = positive_branch[0]
    ax.scatter(
        [gamma["kappa"]],
        [gamma["Omega_real"]],
        s=22,
        facecolor="#6D4C8D",
        edgecolor=WHITE,
        linewidth=0.6,
        zorder=10,
    )
    # The asymmetric two-groove cell does not establish a mirror-symmetry
    # protection rule at Gamma.  Use the neutral label "dark mode" here;
    # the strict Rayleigh-BIC claim is reserved for the red threshold point.
    ax.text(
        0.018,
        gamma["Omega_real"] - 0.035,
        r"$\Gamma$ dark mode",
        color="#6D4C8D",
        fontsize=6.5,
        ha="left",
        va="top",
    )

    strict_id = np.where(values["Omega_imag"] == 0.0)[0]
    assert len(strict_id) == 1
    k_bic = float(values["kappa"][strict_id[0]])
    omega_bic = float(values["Omega_real"][strict_id[0]])
    ax.scatter(
        [-k_bic, k_bic],
        [omega_bic, omega_bic],
        s=25,
        facecolor=RED,
        edgecolor=WHITE,
        linewidth=0.6,
        zorder=10,
    )
    ax.text(0.135, 0.925, "Rayleigh BIC", color=RED, fontsize=6.7, ha="left")
    ax.text(0.335, 0.355, r"$n=0$", color="#77838B", fontsize=6.5, rotation=42, rotation_mode="anchor")
    ax.text(0.305, 0.665, r"$n=-1$", color=ORANGE, fontsize=6.5, rotation=-42, rotation_mode="anchor")
    ax.text(-0.305, 0.665, r"$n=+1$", color=ORANGE, fontsize=6.5, rotation=42, rotation_mode="anchor", ha="right")
    ax.text(
        0.265,
        0.115,
        "guided",
        color="#667179",
        fontsize=6.4,
        ha="left",
        va="center",
    )

    ax.text(
        -0.475,
        0.055,
        "outgoing-pole spectrum",
        color=MUTED,
        fontsize=6.1,
        ha="left",
        va="bottom",
    )
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0.0, 1.0)
    ax.set_xticks([-0.5, -0.25, 0.0, 0.25, 0.5])
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel(r"Bloch wave number, $\kappa=k_xa/(2\pi)$")
    ax.set_ylabel(r"$\mathrm{Re}\,\Omega_p$")
    style_axis(ax)

    scalar = mpl.cm.ScalarMappable(norm=Q_NORM, cmap=Q_CMAP)
    scalar.set_array([])
    cb = plt.colorbar(scalar, cax=cax)
    cb.set_ticks([0, 3, 6, 9, 12])
    cb.set_ticklabels(["0", "3", "6", "9", "≥12"])
    cb.ax.tick_params(labelsize=6.4, width=0.6, length=2.2, pad=1.5)
    cb.outline.set_linewidth(0.6)
    cb.set_label(r"$\log_{10}Q$", fontsize=7.2, labelpad=2.0)


def load_evolution() -> np.ndarray:
    values = np.genfromtxt(EVOLUTION_CSV, delimiter=",", names=True)
    values = values[np.argsort(values["kappa"])]
    endpoint = int(np.argmin(np.abs(values["delta"])))
    assert values["A0_over_C"][endpoint] == 0.0
    assert values["Am1_over_C"][endpoint] == 0.0
    return values


def load_target_bz() -> np.ndarray:
    values = np.genfromtxt(TARGET_BZ_CSV, delimiter=",", names=True)
    values = values[np.argsort(values["kappa"])]
    assert values["kappa"][0] <= 1e-12
    assert abs(values["kappa"][-1] - 0.5) <= 1e-12
    assert np.nanmax(values["sigma_ratio"]) < 1e-9
    return values


def draw_zoomed_pole(ax: plt.Axes, values: np.ndarray, g: dict[str, float]) -> None:
    kappa = values["kappa"]
    omega = values["Omega_real"]
    imag = np.abs(values["Omega_imag"])
    q = np.divide(
        omega,
        2.0 * imag,
        out=np.full_like(omega, np.nan),
        where=imag > 0,
    )
    window = (kappa >= 0.103) & (kappa <= 0.117)
    ax.plot(kappa[window], omega[window], color=NAVY, linewidth=1.45, zorder=4)
    line_k = np.linspace(0.103, 0.117, 301)
    ax.plot(
        line_k,
        1.0 - line_k,
        color=ORANGE,
        linewidth=1.05,
        linestyle=(0, (4.0, 2.3)),
        zorder=2,
    )
    ax.axvline(g["kappa"], color=RED, linewidth=0.65, linestyle=(0, (2.2, 2.2)))
    ax.plot(g["kappa"], g["omega"], "o", ms=4.8, mfc=RED, mec=WHITE, mew=0.55, zorder=8)
    ax.set_xlim(0.103, 0.117)
    ax.set_ylim(0.8825, 0.8975)
    ax.set_xticks([0.104, 0.108, 0.112, 0.116])
    ax.set_yticks([0.885, 0.890, 0.895])
    ax.set_xlabel(r"$\kappa$")
    ax.set_ylabel(r"$\mathrm{Re}\,\Omega_p$")
    style_axis(ax, grid=True)

    axq = ax.twinx()
    root = int(np.argmin(np.abs(values["delta"])))
    for side in (slice(0, root), slice(root + 1, None)):
        valid = window[side] & np.isfinite(q[side])
        axq.semilogy(
            kappa[side][valid],
            q[side][valid],
            color=PLUM,
            linewidth=1.05,
            marker="o",
            markersize=1.9,
            markerfacecolor=WHITE,
            markeredgecolor=PLUM,
            markeredgewidth=0.45,
            zorder=3,
        )
    axq.set_ylim(1e5, 1e15)
    axq.set_ylabel(r"$Q$", color=PLUM, labelpad=2.0)
    axq.tick_params(axis="y", colors=PLUM, pad=2.0, labelsize=7.2)
    axq.spines["right"].set_visible(True)
    axq.spines["right"].set_color(PLUM)
    axq.spines["top"].set_visible(False)
    axq.annotate(
        r"$Q\rightarrow\infty$",
        xy=(g["kappa"], 4.0e14),
        xytext=(g["kappa"] + 0.0013, 1.5e14),
        fontsize=6.5,
        color=PLUM,
        ha="left",
        arrowprops=dict(arrowstyle="->", color=PLUM, linewidth=0.6),
    )
    handles = [
        Line2D([0], [0], color=NAVY, linewidth=1.4, label=r"$\mathrm{Re}\,\Omega_p$"),
        Line2D([0], [0], color=ORANGE, linewidth=1.0, linestyle="--", label=r"$n=-1$ Rayleigh line"),
        Line2D([0], [0], color=PLUM, linewidth=1.0, marker="o", markersize=2.2, markerfacecolor=WHITE, label=r"$Q$"),
    ]
    ax.legend(
        handles=handles,
        loc="lower left",
        fontsize=7.2,
        handlelength=1.8,
        handletextpad=0.45,
        labelspacing=0.28,
        borderaxespad=0.35,
    )


def draw_channel_closure(ax: plt.Axes, values: np.ndarray, g: dict[str, float]) -> None:
    kappa = values["kappa"]
    a0 = values["A0_over_C"]
    am1 = values["Am1_over_C"]
    root = int(np.argmin(np.abs(values["delta"])))
    ax.axvspan(kappa.min(), g["kappa"], color="#E9EDF0", alpha=0.72, zorder=0)
    ax.axvspan(g["kappa"], kappa.max(), color="#FBECDD", alpha=0.58, zorder=0)
    for side in (slice(0, root), slice(root + 1, None)):
        valid = a0[side] > 0
        ax.semilogy(
            kappa[side][valid], a0[side][valid], color=BLUE,
            linewidth=1.25, marker="o", markersize=2.35,
            markerfacecolor=WHITE, markeredgecolor=BLUE, markeredgewidth=0.5,
            zorder=4,
        )
    left = slice(0, root)
    right = slice(root + 1, None)
    ax.semilogy(
        kappa[left], am1[left], color=MUTED, linewidth=1.2,
        linestyle=(0, (3.3, 2.2)), marker="o", markersize=2.25,
        markerfacecolor=WHITE, markeredgecolor=MUTED, markeredgewidth=0.5,
        zorder=4,
    )
    ax.semilogy(
        kappa[right], am1[right], color=ORANGE, linewidth=1.25,
        marker="o", markersize=2.35, markerfacecolor=WHITE,
        markeredgecolor=ORANGE, markeredgewidth=0.5, zorder=4,
    )
    ax.axvline(g["kappa"], color=RED, linewidth=0.65, linestyle=(0, (2.2, 2.2)))
    ax.annotate(
        r"$A_0=A_{-1}=0$",
        xy=(g["kappa"], 0.0),
        xycoords=("data", "axes fraction"),
        xytext=(g["kappa"] + 0.0030, 2.4e-7),
        textcoords="data",
        fontsize=7.2,
        color=RED,
        arrowprops=dict(arrowstyle="->", color=RED, linewidth=0.6),
    )
    ax.text(0.0975, 3.0e-2, r"$n=-1$ closed", color=MUTED, fontsize=6.5, ha="center")
    ax.text(0.1225, 3.0e-2, r"$n=-1$ open", color=ORANGE, fontsize=6.5, ha="center")
    ax.text(0.132, 1.20e-2, r"$|A_{-1}|$", color=ORANGE, fontsize=7.2, ha="right")
    ax.text(0.132, 2.1e-3, r"$|A_0|$", color=BLUE, fontsize=7.2, ha="right")
    ax.set_xlim(0.085, 0.135)
    ax.set_ylim(3e-8, 4e-2)
    ax.set_xticks([0.09, 0.10, 0.11, 0.12, 0.13])
    ax.set_xlabel(r"$\kappa$")
    ax.set_ylabel(r"Floquet coefficient, $|A_n|/\Vert C\Vert_2$")
    style_axis(ax, grid=True)


def load_fields():
    raw = loadmat(FIELD_MAT, squeeze_me=True, struct_as_record=False)
    states = list(np.atleast_1d(raw["states"]))
    x = np.asarray(raw["x"], dtype=float)
    y = np.asarray(raw["y"], dtype=float)
    assert len(states) == 3
    assert float(states[0].sigma_ratio) < 1e-12
    assert float(states[1].sigma_ratio) < 1e-12
    assert float(states[2].raw_residual) < 1e-9
    return x, y, states


def classify_channel(state, order: int, ky: complex) -> str:
    scale = max(abs(complex(state.Omega)) * 2.0 * np.pi, 1.0)
    if abs(ky) <= 2e-7 * scale:
        return "Rayleigh"
    if real_positive(ky) and abs(np.imag(ky)) <= 1e-3 * max(abs(np.real(ky)), 1.0):
        return "propagating"
    return "evanescent" if np.imag(ky) < 0 else "adjacent sheet"


def real_positive(value: complex) -> bool:
    return np.real(value) > 0


def state_channel(state, order: int) -> tuple[complex, complex, complex]:
    suffix = {1: "p1", 0: "0", -1: "m1"}[order]
    return (
        complex(getattr(state, f"A{suffix}_complex")),
        complex(getattr(state, f"kx{suffix}")),
        complex(getattr(state, f"ky{suffix}")),
    )


def add_wavevector_arrow(
    ax: plt.Axes, kx: complex, ky: complex, amplitude: float, status: str
) -> None:
    if amplitude == 0.0:
        ax.text(
            0.5,
            0.48,
            "×",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color=RED,
            fontsize=15,
            fontweight="bold",
        )
        return
    if status != "propagating":
        return
    dx = float(np.real(kx))
    dy = float(np.real(ky))
    norm = max(np.hypot(dx, dy), np.finfo(float).eps)
    dx, dy = 0.34 * dx / norm, 0.34 * dy / norm
    cx, cy = 0.50, 0.50
    ax.add_patch(
        FancyArrowPatch(
            (cx - 0.5 * dx, cy - 0.5 * dy),
            (cx + 0.5 * dx, cy + 0.5 * dy),
            transform=ax.transAxes,
            arrowstyle="-|>",
            mutation_scale=7.0,
            linewidth=1.0,
            linestyle="-",
            color=WHITE,
            zorder=8,
        )
    )


def draw_channel_decomposition(
    axes: list[list[plt.Axes]], cax: plt.Axes, g: dict[str, float]
) -> None:
    del g  # The decomposition is evaluated only in the homogeneous exterior.
    _, _, states = load_fields()
    orders = [1, 0, -1]
    headers = [r"$\Gamma$ dark mode", "Rayleigh BIC", r"$\kappa>\kappa_{\mathrm{BIC}}$"]
    # A tall exterior window makes the channel topology visible without any
    # a-priori open/closed filtering: closed orders decay naturally, whereas
    # propagating orders survive into the far field.
    x = np.linspace(0.0, 5.0, 1001)
    y = np.linspace(0.0, 3.2, 641)
    xx, yy = np.meshgrid(x, y)
    image = None

    for row, order in enumerate(orders):
        for col, state in enumerate(states):
            ax = axes[row][col]
            amplitude, kx, ky = state_channel(state, order)
            magnitude = abs(amplitude)
            cancelled = magnitude < 1e-10
            if cancelled:
                displayed_field = np.zeros((y.size, x.size))
            else:
                channel_field = amplitude * np.exp(-1j * kx * xx - 1j * ky * yy)
                # Normalize only by the pressure magnitude at the surface.
                # This preserves the actual vertical decay/propagation law;
                # it does not equalize a closed channel at large y.
                displayed_field = np.real(channel_field / magnitude)
            image = ax.imshow(
                displayed_field,
                origin="lower",
                extent=[0.0, 5.0, 0.0, 3.2],
                cmap=CHANNEL_REAL_CMAP,
                vmin=-1.0,
                vmax=1.0,
                interpolation="nearest",
                aspect="auto",
                rasterized=True,
            )
            status = classify_channel(state, order, ky)
            add_wavevector_arrow(ax, kx, ky, 0.0 if cancelled else magnitude, status)
            order_text = {1: "+1", 0: "0", -1: "-1"}[order]
            if cancelled and order == 0:
                status = "open · cancelled"
            elif cancelled and order == -1:
                status = "Rayleigh · cancelled"
            ax.set_title(
                rf"$n={order_text}$   " + status,
                pad=2.0,
                fontsize=7.2,
                color=RED if cancelled else INK,
            )
            if row == 0:
                ax.text(
                    0.5,
                    1.23,
                    headers[col],
                    transform=ax.transAxes,
                    ha="center",
                    va="bottom",
                    fontsize=7.4,
                    color=INK,
                    clip_on=False,
                )
            if col == 0:
                ax.set_yticks([0.0, 1.5, 3.0])
            else:
                ax.set_yticks([0.0, 1.5, 3.0])
                ax.tick_params(labelleft=False)
            if row == len(orders) - 1:
                ax.set_xticks([0.0, 2.5, 5.0])
                ax.set_xlabel(r"$x/a$")
            else:
                ax.set_xticks([0.0, 2.5, 5.0])
                ax.tick_params(labelbottom=False)
            ax.set_xlim(0.0, 5.0)
            ax.set_ylim(0.0, 3.2)
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_color(INK)
                spine.set_linewidth(0.55)
            ax.tick_params(length=2.0, width=0.55, pad=1.2)

    axes[1][0].set_ylabel(r"$y/a$", labelpad=2.0)
    assert image is not None
    cb = plt.colorbar(image, cax=cax)
    cb.set_ticks([-1, 0, 1])
    cb.ax.tick_params(labelsize=6.6, width=0.55, length=2.2, pad=1.5)
    cb.outline.set_linewidth(0.55)
    cb.set_label(
        r"$\mathrm{Re}\,p_n/|p_n(y=0)|$", fontsize=7.2, labelpad=2.0
    )


def write_manifest(g: dict[str, float]) -> None:
    out = DATA / "fig2_three_state_fields_180k" / "fig2_source_manifest.txt"
    lines = [
        "Fig. 2 source manifest",
        "",
        "(a) Geometry: StrictRayleighBIC_180kHz_7p10deg_final.mat",
        f"    a={g['a_mm']:.9g} mm, w1={g['w1_mm']:.9g} mm, "
        f"w2={g['w2_mm']:.9g} mm, g={g['gap_mm']:.9g} mm, "
        f"d1={g['d1_mm']:.9g} mm, d2={g['d2_mm']:.9g} mm.",
        "(b) Bands: physical-boundary branches plus the directly recomputed target",
        "    eigenvalue branch from Gamma to the BZ edge, mirrored by reciprocity.",
        "    Curve position is Re(Omega_p); colour is log10(Q), and therefore",
        "    encodes the vanishing imaginary part. The same target branch contains",
        "    a Gamma dark mode and the off-Gamma Rayleigh BIC.",
        "    Adjacent-sheet branches are",
        "    excluded from the main panel.",
        "(c,d) kx_channel_evolution_physical_180k/physical_sheet_channel_evolution.csv",
        "    Every off-BIC row satisfies the all-channel outgoing/decaying operator;",
        "    the exact strict endpoint contains A0=A-1=0 and ImOmega=0.",
        "(e) fig2_three_state_fields_180k/eigenmode_fields_gamma_bic_beyond.mat",
        "    The physical-sheet homogeneous exterior field is decomposed into",
        "    n=+1,0,-1 Floquet components over five periods and 3.2 periods in",
        "    height. Each nonzero component is normalized only by its magnitude at",
        "    y=0, preserving its actual vertical decay or propagation law.",
        "    Exact zeros remain exact and are shown by a red cross.",
        "    Gamma is shown only as the single-open-channel dark mode; no",
        "    mirror-symmetry protection is assumed for the asymmetric cell.",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_figure() -> plt.Figure:
    g = load_geometry()
    values = load_evolution()
    target_bz = load_target_bz()
    write_manifest(g)

    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM), facecolor=WHITE)
    outer = fig.add_gridspec(
        3,
        1,
        height_ratios=[1.05, 0.72, 2.05],
        left=0.052,
        right=0.965,
        bottom=0.065,
        top=0.975,
        hspace=0.34,
    )
    top = outer[0].subgridspec(1, 2, width_ratios=[0.37, 0.63], wspace=0.22)
    right_top = top[0, 1].subgridspec(1, 2, width_ratios=[1.0, 0.035], wspace=0.10)
    middle = outer[1].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.36)
    bottom = outer[2].subgridspec(
        3,
        4,
        width_ratios=[1.0, 1.0, 1.0, 0.035],
        height_ratios=[1.0, 1.0, 1.0],
        wspace=0.12,
        hspace=0.25,
    )

    ax_a = fig.add_subplot(top[0, 0])
    ax_b = fig.add_subplot(right_top[0, 0])
    cax_b = fig.add_subplot(right_top[0, 1])
    ax_c = fig.add_subplot(middle[0, 0])
    ax_d = fig.add_subplot(middle[0, 1])
    field_axes = [
        [fig.add_subplot(bottom[i, j]) for j in range(3)] for i in range(3)
    ]
    cax_e = fig.add_subplot(bottom[:, 3])

    draw_geometry(ax_a, g)
    draw_full_band(ax_b, cax_b, target_bz)
    draw_zoomed_pole(ax_c, values, g)
    draw_channel_closure(ax_d, values, g)
    draw_channel_decomposition(field_axes, cax_e, g)

    panel_label(ax_a, "a", x=-0.015, y=1.015)
    panel_label(ax_b, "b", x=-0.105, y=1.035)
    panel_label(ax_c, "c", x=-0.105, y=1.045)
    panel_label(ax_d, "d", x=-0.105, y=1.045)
    panel_label(field_axes[0][0], "e", x=-0.115, y=1.25)
    return fig


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig = make_figure()
    base = FIG_DIR / "fig2_rayleigh_bic_eigenproblem"
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.025)
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.025)
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight", pad_inches=0.025)
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pad_inches=0.025)
    plt.close(fig)


if __name__ == "__main__":
    main()
