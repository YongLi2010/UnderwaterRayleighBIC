#!/usr/bin/env python3
"""Corrected full complex-band panel for the 180-kHz Rayleigh BIC.

The MATLAB calculations remain the numerical source.  This script only
classifies and renders those poles: states below the n=0 sound line are shown
as guided/bound and are deliberately excluded from the Q colour scale.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "Ni2019_MATLAB" / "results" / "full_complex_band_180k"
FIG_DIR = ROOT / "arxiv_theory_paper" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

MM = 1.0 / 25.4
WIDTH_MM = 178.0
HEIGHT_MM = 104.0

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7.5,
        "axes.labelsize": 8.5,
        "axes.linewidth": 0.75,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "xtick.labelsize": 7.2,
        "ytick.labelsize": 7.2,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    }
)

INK = "#20252A"
NEUTRAL = "#7A858E"
NEUTRAL_LIGHT = "#F1F3F4"
SOUND_LINE = "#87929B"
RAYLEIGH = "#D97716"
ACCENT = "#D43F34"

Q_MAX = 12.0
Q_CMAP = LinearSegmentedColormap.from_list(
    "quiet_to_radiant_q",
    ["#C9D0D4", "#80AFC1", "#278F98", "#F0B23A", "#DC4638"],
)
Q_NORM = Normalize(vmin=0.0, vmax=Q_MAX, clip=True)

BRANCH_NAMES = {
    1: "broad_resonance_I",
    2: "narrow_resonance",
    4: "guided_below_sound_line",
    9: "broad_resonance_II",
    10: "broad_resonance_I_adjacent_sheet",
    19: "broad_resonance_II_adjacent_sheet",
    100: "target_rayleigh_sheet",
}


def read_numeric_csv(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        return [
            {key: float(value) for key, value in row.items() if value != ""}
            for row in reader
        ]


def load_source_data() -> list[dict[str, Any]]:
    raw: list[dict[str, Any]] = read_numeric_csv(
        DATA_DIR / "full_complex_band_points.csv"
    )
    for row in raw:
        row["track_id"] = int(row["track_id"])
        row["source"] = "independent pole census / target-sheet continuation"
        row["sheet_role"] = (
            "target Rayleigh sheet" if row["track_id"] == 100 else "physical sheet"
        )
        row["full_operator_sigma_ratio"] = np.nan
        row["full_operator_raw_residual"] = np.nan
        row["continuation_sigma_ratio"] = np.nan
        row["continuation_residual"] = np.nan
        row["mode_overlap"] = np.nan

    near: list[dict[str, Any]] = []
    for source_row in read_numeric_csv(DATA_DIR / "low_guided_branch_near_gamma.csv"):
        row = {
            key: source_row[key]
            for key in ("kappa", "Omega_real", "Omega_imag", "Q")
        }
        row.update(
            {
                "track_id": 4,
                "source": "direct near-Gamma continuation",
                "sheet_role": "below-continuum guided state",
                "full_operator_sigma_ratio": np.nan,
                "full_operator_raw_residual": np.nan,
                "continuation_sigma_ratio": np.nan,
                "continuation_residual": np.nan,
                "mode_overlap": np.nan,
            }
        )
        near.append(row)
    mirrored = [dict(row, kappa=-float(row["kappa"])) for row in near]

    strict_row = read_numeric_csv(DATA_DIR / "strict_bic_endpoint.csv")[0]
    strict: dict[str, Any] = dict(strict_row)
    strict.update(
        {
            "track_id": 100,
            "Q": np.inf,
            "source": "strict full-operator root",
            "sheet_role": "target Rayleigh sheet",
            "continuation_sigma_ratio": np.nan,
            "continuation_residual": np.nan,
            "mode_overlap": np.nan,
        }
    )
    reciprocal = dict(
        strict,
        kappa=-float(strict["kappa"]),
        source="strict full-operator root (reciprocal partner)",
    )

    k_bic = float(strict["kappa"])
    raw = [
        row
        for row in raw
        if not (
            row["track_id"] == 100
            and abs(abs(row["kappa"]) - k_bic) < 5e-10
            and abs(row["Omega_imag"]) < 1e-18
        )
    ]

    adjacent: list[dict[str, Any]] = []
    sheet_specs = [
        ("dark_branch_n0_fixed_sheet.csv", 10, "fixed n=0 sheet continuation"),
        ("dark_branch_nm1_fixed_sheet.csv", 19, "fixed n=-1 sheet continuation"),
    ]
    for filename, track_id, label in sheet_specs:
        for source_row in read_numeric_csv(DATA_DIR / filename):
            row = {
                "track_id": track_id,
                "kappa": source_row["kappa"],
                "Omega_real": source_row["Omega_real"],
                "Omega_imag": source_row["Omega_imag"],
                "Q": source_row["Q"],
                "source": label,
                "sheet_role": "adjacent-sheet continuation",
                "full_operator_sigma_ratio": np.nan,
                "full_operator_raw_residual": np.nan,
                "continuation_sigma_ratio": source_row["sigma_ratio"],
                "continuation_residual": source_row["residual"],
                "mode_overlap": source_row["mode_overlap"],
            }
            adjacent.append(row)
            if abs(row["kappa"]) > 1e-14:
                mirrored_row = dict(row, kappa=-float(row["kappa"]))
                mirrored_row["source"] = label + " (reciprocal partner)"
                adjacent.append(mirrored_row)

    deduplicated: dict[tuple[int, float], dict[str, Any]] = {}
    for row in [*raw, *near, *mirrored, *adjacent, strict, reciprocal]:
        key = (int(row["track_id"]), round(float(row["kappa"]), 12))
        deduplicated.setdefault(key, row)
    data = sorted(
        deduplicated.values(), key=lambda row: (row["track_id"], row["kappa"])
    )

    for row in data:
        row["branch"] = BRANCH_NAMES.get(int(row["track_id"]), "")
        row["open_channel_count"] = int(
            sum(
                row["Omega_real"] > abs(row["kappa"] + order)
                for order in range(-3, 4)
            )
        )
        if row["track_id"] in (10, 19):
            row["state_type"] = "adjacent-sheet pole"
        elif row["open_channel_count"] == 0:
            row["state_type"] = "guided/bound"
        else:
            row["state_type"] = "radiative pole"
        is_bic = (
            row["track_id"] == 100
            and abs(row["Omega_imag"]) < 1e-18
            and abs(row["Omega_real"] - (1.0 - abs(row["kappa"]))) < 5e-10
        )
        if is_bic:
            row["state_type"] = "strict Rayleigh BIC"
        row["log10_Q_display"] = (
            float(
                np.clip(
                    np.log10(np.clip(float(row["Q"]), 1.0, np.inf)),
                    0.0,
                    Q_MAX,
                )
            )
            if row["open_channel_count"] > 0 or row["track_id"] in (10, 19)
            else np.nan
        )
        row["displayed"] = bool(
            0.0 <= row["Omega_real"] <= 1.0 and row["branch"]
        )
    return data


def split_contiguous(
    branch: list[dict[str, Any]],
) -> list[list[dict[str, Any]]]:
    branch = sorted(branch, key=lambda row: row["kappa"])
    if len(branch) < 2:
        return [branch]
    dx = np.diff([row["kappa"] for row in branch])
    positive = dx[dx > 1e-8]
    nominal = np.median(positive) if len(positive) else 0.01
    cut = np.where(dx > 1.6 * nominal)[0] + 1
    indices = [0, *cut.tolist(), len(branch)]
    return [branch[start:stop] for start, stop in zip(indices[:-1], indices[1:])]


def add_q_coloured_curve(
    ax: plt.Axes,
    branch: list[dict[str, Any]],
    linewidth: float,
    linestyle: str = "solid",
    alpha: float = 1.0,
) -> None:
    for part in split_contiguous(branch):
        if len(part) < 2:
            continue
        x = np.asarray([row["kappa"] for row in part], dtype=float)
        y = np.asarray([row["Omega_real"] for row in part], dtype=float)
        q = np.asarray([row["log10_Q_display"] for row in part], dtype=float)
        points = np.column_stack([x, y])
        segments = np.stack([points[:-1], points[1:]], axis=1)
        values = 0.5 * (q[:-1] + q[1:])
        collection = LineCollection(
            segments,
            cmap=Q_CMAP,
            norm=Q_NORM,
            linewidth=linewidth,
            linestyles=linestyle,
            alpha=alpha,
            capstyle="round",
            joinstyle="round",
            zorder=4,
        )
        collection.set_array(values)
        ax.add_collection(collection)


def add_plain_curve(
    ax: plt.Axes,
    branch: list[dict[str, Any]],
    linewidth: float,
    color: str,
    alpha: float,
) -> None:
    for part in split_contiguous(branch):
        if len(part) < 2:
            continue
        ax.plot(
            [row["kappa"] for row in part],
            [row["Omega_real"] for row in part],
            color=color,
            linewidth=linewidth,
            linestyle=(0, (3.2, 2.2)),
            dash_capstyle="butt",
            alpha=alpha,
            zorder=3.7,
        )


def make_figure(data: list[dict[str, Any]]) -> plt.Figure:
    fig, ax = plt.subplots(
        figsize=(WIDTH_MM * MM, HEIGHT_MM * MM), constrained_layout=False
    )
    fig.subplots_adjust(left=0.105, right=0.855, bottom=0.16, top=0.965)

    kappa = np.linspace(-0.5, 0.5, 1001)
    ax.fill_between(
        kappa,
        0.0,
        np.abs(kappa),
        color=NEUTRAL_LIGHT,
        linewidth=0,
        zorder=0,
    )
    ax.plot(
        kappa,
        np.abs(kappa),
        color=SOUND_LINE,
        linewidth=0.9,
        linestyle=(0, (4.0, 3.0)),
        zorder=2,
    )
    ax.plot(
        kappa,
        1.0 - np.abs(kappa),
        color=RAYLEIGH,
        linewidth=1.05,
        linestyle=(0, (5.0, 3.2)),
        zorder=2,
    )

    visible = [row for row in data if row["displayed"]]
    guided = [row for row in visible if row["track_id"] == 4]
    for part in split_contiguous(guided):
        ax.plot(
            [row["kappa"] for row in part],
            [row["Omega_real"] for row in part],
            color=NEUTRAL,
            linewidth=1.35,
            solid_capstyle="round",
            zorder=3,
        )

    for track_id in (1, 2, 9, 100):
        branch = [row for row in visible if row["track_id"] == track_id]
        add_q_coloured_curve(
            ax, branch, linewidth=1.65 if track_id == 100 else 1.15
        )

    for track_id in (10, 19):
        branch = [row for row in visible if row["track_id"] == track_id]
        add_plain_curve(
            ax, branch, linewidth=0.95, color="#43515B", alpha=0.78
        )

    bic = [row for row in visible if row["state_type"] == "strict Rayleigh BIC"]
    ax.scatter(
        [row["kappa"] for row in bic],
        [row["Omega_real"] for row in bic],
        s=34,
        facecolor=ACCENT,
        edgecolor="white",
        linewidth=0.8,
        zorder=8,
        clip_on=False,
    )

    ax.text(
        0.137,
        0.925,
        "Rayleigh BIC",
        color=ACCENT,
        fontsize=7.4,
        ha="left",
        va="bottom",
    )
    ax.text(
        0.30,
        0.155,
        "guided band\n(no open channel)",
        color="#626D75",
        fontsize=7.2,
        ha="left",
        va="center",
        linespacing=1.05,
    )
    sheet_handles = [
        Line2D([0], [0], color="#43515B", linewidth=1.1, linestyle="solid"),
        Line2D([0], [0], color="#43515B", linewidth=1.0, linestyle="dashed"),
    ]
    ax.legend(
        sheet_handles,
        ["physical sheet", "adjacent sheet"],
        loc="lower left",
        bbox_to_anchor=(0.012, 0.018),
        borderaxespad=0,
        handlelength=2.4,
        handletextpad=0.6,
        labelspacing=0.35,
        fontsize=7.2,
    )
    ax.text(
        0.365,
        0.440,
        r"$n=0$",
        color=SOUND_LINE,
        fontsize=7.2,
        rotation=40,
        rotation_mode="anchor",
        ha="left",
        va="bottom",
    )
    ax.text(
        0.315,
        0.615,
        r"$n=-1$",
        color=RAYLEIGH,
        fontsize=7.2,
        rotation=-40,
        rotation_mode="anchor",
        ha="left",
        va="bottom",
    )
    ax.text(
        -0.315,
        0.615,
        r"$n=+1$",
        color=RAYLEIGH,
        fontsize=7.2,
        rotation=40,
        rotation_mode="anchor",
        ha="right",
        va="bottom",
    )

    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0.0, 1.0)
    ax.set_xticks([-0.5, -0.25, 0.0, 0.25, 0.5])
    ax.set_yticks(np.linspace(0.0, 1.0, 6))
    ax.set_xlabel(r"Normalized Bloch wave number, $\kappa=k_xa/(2\pi)$")
    ax.set_ylabel(r"Normalized eigenfrequency, $\mathrm{Re}\,\Omega$")
    ax.tick_params(colors=INK)
    ax.spines["left"].set_color(INK)
    ax.spines["bottom"].set_color(INK)
    ax.set_facecolor("white")

    sm = mpl.cm.ScalarMappable(norm=Q_NORM, cmap=Q_CMAP)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.038, pad=0.025, aspect=24)
    cbar.set_ticks([0, 3, 6, 9, 12])
    cbar.set_ticklabels(["0", "3", "6", "9", "≥12"])
    cbar.ax.tick_params(labelsize=7.2, width=0.65, length=2.5)
    cbar.outline.set_linewidth(0.65)
    cbar.set_label(r"$\log_{10}Q$", fontsize=7.5, labelpad=4)

    return fig


def write_source_data(data: list[dict[str, Any]]) -> None:
    fields = [
        "track_id",
        "branch",
        "kappa",
        "Omega_real",
        "Omega_imag",
        "Q",
        "open_channel_count",
        "state_type",
        "log10_Q_display",
        "displayed",
        "source",
        "full_operator_sigma_ratio",
        "full_operator_raw_residual",
        "continuation_sigma_ratio",
        "continuation_residual",
        "mode_overlap",
        "sheet_role",
    ]
    with (DATA_DIR / "fig2b_full_complex_band_source_data.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data)


def write_qa_notes(data: list[dict[str, Any]]) -> None:
    shown = [row for row in data if row["displayed"]]
    guided = [row for row in shown if row["state_type"] == "guided/bound"]
    adjacent = [
        row for row in shown if row["state_type"] == "adjacent-sheet pole"
    ]
    lines = [
        "# Fig. 2b corrected full-band panel: QA notes",
        "",
        "- Core conclusion: an outgoing pole crosses the Rayleigh line while its radiative linewidth vanishes; the low-frequency grey branch is a separate below-continuum guided state.",
        "- Archetype: quantitative discovery panel.",
        f"- Displayed points: {len(shown)}; source rows retained before the 0–1 window filter: {len(data)}.",
        f"- Guided/bound points excluded from Q colour mapping: {len(guided)}.",
        f"- Fixed-sheet continuation points shown with dashed lines: {len(adjacent)}.",
        "- Display window: 0 ≤ Re Omega ≤ 1, as requested.",
        "- Branches are split whenever the kappa gap exceeds 1.6 times the local sampling step; no line is drawn across a missing Riemann-sheet segment.",
        "- The two broad low-Q branches are distinct. Dashed segments are q-parameterized continuations on the adjacent sheet and are never joined to each other.",
        "- Q colours are clipped to 0 ≤ log10(Q) ≤ 12. Low-Q poles are deliberately quiet, whereas high-Q states are vivid; exact BICs are identified independently by markers.",
        "- The near-Gamma guided points were recalculated in MATLAB and are identified in the source column.",
        "- No interpolation, smoothing, or simulated data are used.",
    ]
    (DATA_DIR / "fig2b_full_complex_band_QA.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    data = load_source_data()
    write_source_data(data)
    write_qa_notes(data)
    fig = make_figure(data)
    base = FIG_DIR / "fig2b_full_complex_band"
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
