#!/usr/bin/env python3
"""Fig. 4 | Experimental observability of the 180-kHz Rayleigh BIC.

Every quantitative value is theory.  Continuous maps and lines show the full
theoretical prediction.  Open symbols and the coarse field sampling are exact
subsamples of the same theory and reserve the future experimental data layer;
they must not be described as measurements until replaced by measured data.
"""
from pathlib import Path
import csv

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Arc, Circle, FancyArrowPatch, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DATA = (ROOT / "Ni2019_MATLAB" / "results" /
        "fig4_experimental_observables_180k")
OUT = ROOT / "arxiv_theory_paper" / "figures"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.4,
    "axes.titlesize": 8.1,
    "axes.titleweight": "semibold",
    "axes.labelsize": 7.4,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.65,
    "ytick.major.width": 0.65,
    "xtick.major.size": 2.6,
    "ytick.major.size": 2.6,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

INK = "#24333B"
NAVY = "#214A65"
TEAL = "#138C8C"
GOLD = "#D8891C"
PURPLE = "#7868A6"
RED = "#D94B3D"
GREY = "#77858B"
LIGHT_GREY = "#D8E0E3"
WATER = "#EAF4F5"
PALE_TEAL = "#D8EEEE"
PALE_GOLD = "#F5E8D0"

PHASE_CMAP = LinearSegmentedColormap.from_list(
    "phase", ["#294C73", "#93BFD1", "#F7F7F4", "#F0B57B", "#A63C32"])
TEAL_CMAP = LinearSegmentedColormap.from_list(
    "teal_response", ["#F7F9F8", "#D8EEEE", "#64B7B5", "#137F82", "#173A55"])
GOLD_CMAP = LinearSegmentedColormap.from_list(
    "gold_response", ["#FAF9F5", "#F5E8D0", "#EABF74", "#D8891C", "#704216"])
FIELD_CMAP = LinearSegmentedColormap.from_list(
    "field", ["#234E72", "#8FC0D2", "#F7F7F4", "#E9A06F", "#A83B33"])


def read_csv(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return {key: np.array([float(row[key]) for row in rows])
            for key in rows[0]}


def panel_label(ax, letter, x=-0.075, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9.3, fontweight="bold", color="#111111", clip_on=False)


def clean_axis(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def style_axis(ax):
    ax.tick_params(labelsize=7.2, pad=1.8)
    ax.spines["left"].set_color(INK)
    ax.spines["bottom"].set_color(INK)


def draw_experiment(ax):
    clean_axis(ax); panel_label(ax, "a", x=-0.03, y=1.01)
    ax.set_title("Water-tank field mapping", pad=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # Subtle tank boundary; white-dominant composition keeps this a diagram,
    # not an illustrative water scene.
    ax.add_patch(Rectangle((0.03, 0.08), 0.94, 0.80,
                           facecolor=WATER, edgecolor="#B8D5DA", lw=0.8))

    # Phase-controlled source array.
    for j in range(7):
        x = 0.11 + 0.041*j
        ax.add_patch(Rectangle((x, 0.72), 0.030, 0.075,
                               angle=-12, facecolor=NAVY,
                               edgecolor="white", lw=0.45))
    ax.text(0.19, 0.82, "phased array", ha="center", color=NAVY)

    # Incident beam and angle marker.
    ax.add_patch(FancyArrowPatch((0.24, 0.69), (0.48, 0.27),
                                 arrowstyle="-|>", mutation_scale=9,
                                 lw=1.15, color=GOLD))
    ax.add_patch(Arc((0.48, 0.27), 0.18, 0.18, theta1=77, theta2=90,
                     color=GOLD, lw=0.8))
    ax.text(0.42, 0.39, r"$\theta$", color=GOLD)

    # Finite reflecting sample. Grooves remain unfilled, matching the actual
    # water-filled recesses rather than being depicted as solid inclusions.
    y0 = 0.24
    ax.plot([0.25, 0.94], [y0, y0], color=INK, lw=1.1)
    for m in range(8):
        xl = 0.29 + 0.078*m
        ax.plot([xl, xl, xl+0.030, xl+0.030],
                [y0, y0-0.105, y0-0.105, y0], color=INK, lw=0.65)
        xs = xl+0.045
        ax.plot([xs, xs, xs+0.014, xs+0.014],
                [y0, y0-0.045, y0-0.045, y0], color=INK, lw=0.65)
    ax.text(0.61, 0.045, "finite metagrating", ha="center", color=INK)

    # Complex-pressure scan window and reference hydrophone.
    ax.add_patch(Rectangle((0.48, 0.34), 0.39, 0.27, fill=False,
                           edgecolor=TEAL, lw=1.0, ls=(0, (3, 2))))
    for xx in np.linspace(0.51, 0.84, 7):
        for yy in np.linspace(0.38, 0.57, 4):
            ax.add_patch(Circle((xx, yy), 0.005, facecolor=TEAL,
                                edgecolor="none", alpha=0.48))
    ax.add_patch(Circle((0.73, 0.50), 0.018, facecolor="white",
                        edgecolor=TEAL, lw=1.05))
    ax.plot([0.73, 0.73], [0.52, 0.62], color=TEAL, lw=0.9)
    ax.text(0.675, 0.65, "complex scan", ha="center", color=TEAL)
    ax.add_patch(Circle((0.87, 0.72), 0.015, facecolor=PURPLE,
                        edgecolor="white", lw=0.45))
    ax.text(0.82, 0.78, "reference", ha="center", color=PURPLE,
            fontsize=7.2)

    # Grazing harmonic is indicated as a field component, not a far-field
    # power arrow.
    ax.add_patch(FancyArrowPatch((0.72, 0.305), (0.31, 0.305),
                                 arrowstyle="-|>", mutation_scale=8,
                                 lw=1.0, color=RED))
    ax.text(0.67, 0.335, r"grazing $n=-1$", color=RED, ha="center")


def assemble_response():
    raw = read_csv(DATA / "adaptive_scattering_response.csv")
    kappa = np.unique(raw["kappa"])
    u = np.unique(raw["normalized_detuning"])
    nk, nu = len(kappa), len(u)
    shape = (nk, nu)
    # MATLAB writes the nK x nU arrays column-by-column.
    def mat(name):
        return raw[name].reshape(shape, order="F")
    theta = mat("theta_deg")[:, 0]
    fp = mat("pole_frequency_hz")[:, 0]
    gamma = mat("linewidth_hz")[:, 0]
    detuning = mat("frequency_hz") - fp[:, None]
    a0 = mat("A0_real") + 1j*mat("A0_imag")
    am = mat("Am1_real") + 1j*mat("Am1_imag")

    # Local complex background from the two remote normalized-detuning ends.
    bg0 = 0.5*(a0[:, 0]+a0[:, -1])
    bgm = 0.5*(am[:, 0]+am[:, -1])
    delta0 = a0-bg0[:, None]
    deltam = am-bgm[:, None]

    y = np.linspace(-4.0, 4.0, 481)
    phase = np.zeros((nk, len(y)))
    response0 = np.zeros_like(phase)
    responsem = np.zeros_like(phase)
    for i in range(nk):
        valid = np.isfinite(a0[i].real) & np.isfinite(am[i].real)
        if valid.sum() < 3 or gamma[i] == 0:
            continue
        x = detuning[i, valid]
        ratio = a0[i, valid]/bg0[i]
        rr = np.interp(y, x, ratio.real, left=1.0, right=1.0)
        ri = np.interp(y, x, ratio.imag, left=0.0, right=0.0)
        phase[i] = np.angle(rr+1j*ri)/np.pi
        response0[i] = np.interp(y, x, abs(delta0[i, valid]),
                                 left=0.0, right=0.0)
        responsem[i] = np.interp(y, x, abs(deltam[i, valid]),
                                 left=0.0, right=0.0)
    return theta, fp, gamma, y, phase, response0, responsem


def heatmap(ax, theta, y, z, cmap, vmin, vmax, title, letter=None):
    image = ax.pcolormesh(theta, y, z.T, shading="auto", cmap=cmap,
                         vmin=vmin, vmax=vmax, rasterized=True)
    ax.set_title(title, pad=4)
    ax.set_xlim(theta.min(), theta.max()); ax.set_ylim(y.min(), y.max())
    ax.set_xlabel(r"Incidence angle $\theta$ (deg)")
    ax.set_ylabel(r"Detuning, $f-\mathrm{Re}\,f_p$ (Hz)")
    style_axis(ax)
    if letter:
        panel_label(ax, letter)
    return image


def draw_linewidth_map(ax, theta, gamma, y, phase):
    im = heatmap(ax, theta, y, phase, PHASE_CMAP, -1, 1,
                 "Collapsing scattering line", "b")
    half = 0.5*gamma
    ax.plot(theta, half, color="white", lw=0.9, zorder=4)
    ax.plot(theta, -half, color="white", lw=0.9, zorder=4)
    ib = int(np.argmin(gamma))
    ax.scatter([theta[ib]], [0], s=27, c=RED, edgecolors="white",
               linewidths=0.55, zorder=6)
    ax.text(theta[ib]+0.075, 3.25, "Rayleigh BIC", color=RED,
            ha="left", va="top", fontsize=7.2)
    # Discrete exact-theory samples reserve future experimental markers.
    ids = np.linspace(5, len(theta)-6, 11).round().astype(int)
    ax.scatter(theta[ids], half[ids], s=15, facecolors="white",
               edgecolors=INK, linewidths=0.55, zorder=5)
    cb = plt.colorbar(im, ax=ax, fraction=0.044, pad=0.025)
    cb.set_label(r"Specular phase $\Delta\phi_0/\pi$")
    cb.ax.tick_params(labelsize=7.2, length=2)


def draw_channels(fig, spec, theta, y, r0, rm):
    outer = fig.add_subplot(spec); clean_axis(outer); panel_label(outer, "c")
    outer.set_title("Channel-resolved resonant response", pad=4)
    sub = spec.subgridspec(2, 1, hspace=0.13)
    ax0 = fig.add_subplot(sub[0]); axm = fig.add_subplot(sub[1], sharex=ax0)

    def log_norm(z):
        scale = np.nanpercentile(z[z > 0], 99.5)
        return np.log10(np.maximum(z/scale, 1e-4))
    z0, zm = log_norm(r0), log_norm(rm)
    im0 = ax0.pcolormesh(theta, y, z0.T, shading="auto", cmap=TEAL_CMAP,
                         vmin=-4, vmax=0, rasterized=True)
    axm.pcolormesh(theta, y, zm.T, shading="auto", cmap=GOLD_CMAP,
                   vmin=-4, vmax=0, rasterized=True)
    for ax in (ax0, axm):
        ax.set_xlim(theta.min(), theta.max()); ax.set_ylim(-4, 4)
        style_axis(ax)
    ax0.tick_params(labelbottom=False)
    ax0.set_ylabel("Detuning (Hz)")
    axm.set_ylabel("Detuning (Hz)")
    axm.set_xlabel(r"Incidence angle $\theta$ (deg)")
    ax0.text(0.02, 0.82, r"$|\Delta r_0|$", transform=ax0.transAxes,
             color=TEAL, fontweight="semibold")
    axm.text(0.02, 0.82, r"$|\Delta r_{-1}|$", transform=axm.transAxes,
             color=GOLD, fontweight="semibold")
    cb = fig.colorbar(im0, ax=[ax0, axm], fraction=0.035, pad=0.018)
    cb.set_label(r"$\log_{10}$ normalized contrast")
    cb.set_ticks([-4, -2, 0]); cb.ax.tick_params(labelsize=7.2, length=2)


def draw_pole_and_residues(fig, spec, summary):
    outer = fig.add_subplot(spec); clean_axis(outer); panel_label(outer, "d")
    outer.set_title("Pole linewidth and channel residues", pad=4)
    sub = spec.subgridspec(2, 1, hspace=0.13)
    ag = fig.add_subplot(sub[0]); ar = fig.add_subplot(sub[1], sharex=ag)
    order = np.argsort(summary["theta_deg"])
    theta = summary["theta_deg"][order]
    gamma = summary["linewidth_hz"][order]
    rad0 = summary["radiation_A0"][order]
    radm = summary["radiation_Am1"][order]
    ib = int(np.argmin(abs(theta-7.099662913)))
    floor_g = 1e-7
    gp = np.maximum(gamma, floor_g)
    ag.semilogy(theta, gp, color=PURPLE, lw=1.45)
    ids = np.linspace(5, len(theta)-6, 11).round().astype(int)
    ag.scatter(theta[ids], gp[ids], s=16, facecolors="white",
               edgecolors=PURPLE, linewidths=0.65, zorder=5)
    ag.scatter(theta[ib], floor_g*1.25, marker="v", s=23, c=RED,
               edgecolors="white", linewidths=0.4, zorder=6)
    ag.set_ylim(floor_g, 10)
    ag.set_ylabel("")
    ag.text(0.018, 0.93, r"$\Gamma$ (Hz)", transform=ag.transAxes,
            ha="left", va="top", color=PURPLE, fontsize=7.4,
            fontweight="semibold")
    ag.text(0.018, 0.10, r"$\circ$ sampled theory", transform=ag.transAxes,
            ha="left", va="bottom", color=INK, fontsize=7.2)
    ag.tick_params(labelbottom=False)

    scale = max(rad0.max(), radm.max())
    r0 = np.maximum(rad0/scale, 1e-15)
    rm = np.maximum(radm/scale, 1e-15)
    ar.semilogy(theta, r0, color=TEAL, lw=1.4)
    ar.semilogy(theta, rm, color=GOLD, lw=1.4, ls=(0, (2.2, 1.5)))
    ar.scatter(theta[ids], r0[ids], s=15, facecolors="white",
               edgecolors=TEAL, linewidths=0.6, zorder=5)
    ar.scatter(theta[ids], rm[ids], s=15, facecolors="white",
               edgecolors=GOLD, linewidths=0.6, zorder=5)
    ar.scatter(theta[ib], max(r0[ib], 1e-15), s=25, c=RED,
               edgecolors="white", linewidths=0.4, zorder=6)
    ar.set_ylim(1e-15, 2)
    ar.set_ylabel("")
    ar.text(0.018, 0.08, r"$|A_n|/\max|A_n|$", transform=ar.transAxes,
            ha="left", va="bottom", color=INK, fontsize=7.4)
    ar.text(0.018, 0.78, r"$n=0$", transform=ar.transAxes,
            color=TEAL, fontsize=7.2)
    ar.text(0.18, 0.78, r"$n=-1$", transform=ar.transAxes,
            color=GOLD, fontsize=7.2)
    ar.set_xlabel(r"Incidence angle $\theta$ (deg)")
    for ax in (ag, ar):
        ax.axvline(theta[ib], color=RED, lw=0.75, ls=(0, (2.4, 2.0)))
        ax.set_xlim(theta.min(), theta.max()); style_axis(ax)

def draw_field_panel(fig, spec):
    outer = fig.add_subplot(spec); clean_axis(outer); panel_label(outer, "e")
    outer.text(0.5, 1.085, "Field-to-Floquet reconstruction",
               transform=outer.transAxes, ha="center", va="bottom",
               fontsize=8.1, fontweight="semibold", color="#111111")
    sub = spec.subgridspec(1, 3, width_ratios=[1.25, 1.25, 0.78], wspace=0.25)
    at = fig.add_subplot(sub[0]); ass = fig.add_subplot(sub[1], sharex=at, sharey=at)
    ak = fig.add_subplot(sub[2])

    x = np.loadtxt(DATA/"field_x_over_a.csv", delimiter=",")
    y = np.loadtxt(DATA/"field_y_over_a.csv", delimiter=",")
    pr = np.loadtxt(DATA/"field_total_real.csv", delimiter=",")
    pi = np.loadtxt(DATA/"field_total_imag.csv", delimiter=",")
    sr = np.loadtxt(DATA/"field_scattered_real.csv", delimiter=",")
    si = np.loadtxt(DATA/"field_scattered_imag.csv", delimiter=",")
    meta = read_csv(DATA/"field_metadata.csv")

    im = at.imshow(pr, origin="lower", aspect="auto",
                   extent=[x.min(), x.max(), y.min(), y.max()],
                   cmap=FIELD_CMAP, vmin=-1, vmax=1, interpolation="bilinear")
    at.set_title("Theory", pad=3, fontsize=7.4)
    at.set_xlabel(r"$x/a$"); at.set_ylabel(r"$y/a$")

    # Coarse exact subsample of the same theory, explicitly not measurement.
    ix = np.arange(0, len(x), 20); iy = np.arange(0, len(y), 16)
    xx, yy = np.meshgrid(x[ix], y[iy])
    ass.scatter(xx.ravel(), yy.ravel(), c=pr[np.ix_(iy, ix)].ravel(),
                cmap=FIELD_CMAP, vmin=-1, vmax=1, s=6.5, marker="s",
                linewidths=0, rasterized=True)
    ass.set_title("Sampled theory", pad=3, fontsize=7.4)
    ass.set_xlabel(r"$x/a$")
    ass.tick_params(labelleft=False)
    for ax in (at, ass):
        ax.set_xlim(-5, 5); ax.set_ylim(0, 3.2); style_axis(ax)
        ax.axhline(0, color=INK, lw=0.7)

    # Quantitative 1D spatial spectrum of the scattered field at y/a ~ 1.5.
    psc = sr+1j*si
    jy = int(np.argmin(abs(y-1.5)))
    window = np.hanning(len(x))
    spectrum = np.fft.fftshift(np.fft.fft(psc[jy]*window))
    dx = x[1]-x[0]
    kx_a = -2*np.pi*np.fft.fftshift(np.fft.fftfreq(len(x), d=dx))
    omega = meta["frequency_hz"][0]*7.416665789e-3/1500.0
    k0_a = 2*np.pi*omega
    kval = kx_a/k0_a
    amp = abs(spectrum); amp /= max(amp)
    mask = (kval > -1.5) & (kval < 1.5)
    order = np.argsort(kval[mask])
    ak.plot(kval[mask][order], amp[mask][order], color=NAVY, lw=1.2)
    theta = np.deg2rad(meta["theta_deg"][0])
    k0 = np.sin(theta); km = k0-1/omega
    for kval0, color, label in [(k0, TEAL, r"$0$"),
                                (km, GOLD, r"$-1$")]:
        ak.axvline(kval0, color=color, lw=0.85, ls=(0, (2.2, 1.5)))
        ak.text(kval0, 0.96, label, color=color, ha="center", va="top")
    ak.set_title("Spatial spectrum", pad=3, fontsize=7.4)
    ak.set_xlabel(r"$k_x/k_0$"); ak.set_ylabel("Amplitude")
    ak.set_xlim(-1.5, 1.5); ak.set_ylim(0, 1.02); style_axis(ak)

    cax = inset_axes(at, width="42%", height="4.0%", loc="upper left",
                     borderpad=0.9)
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_ticks([-1, 0, 1]); cb.ax.tick_params(labelsize=7.2, length=1.7,
                                                pad=0.8)
    cb.set_label(r"$\operatorname{Re}p/\max|p|$", fontsize=7.4, labelpad=0.7)


def main():
    theta, fp, gamma, y, phase, r0, rm = assemble_response()
    summary = read_csv(DATA/"pole_and_radiation_summary.csv")

    assert len(theta) == 91
    assert abs(theta[np.argmin(gamma)]-7.099662913) < 1e-5
    assert np.nanmax(gamma) < 5.0

    fig = plt.figure(figsize=(7.15, 6.35), facecolor="white")
    gs = fig.add_gridspec(3, 12, left=0.055, right=0.985,
                          bottom=0.075, top=0.965,
                          hspace=0.46, wspace=0.55,
                          height_ratios=[0.92, 1.02, 1.05])
    axa = fig.add_subplot(gs[0, :4]); draw_experiment(axa)
    axb = fig.add_subplot(gs[0, 5:]); draw_linewidth_map(axb, theta, gamma, y, phase)
    draw_channels(fig, gs[1, :6], theta, y, r0, rm)
    draw_pole_and_residues(fig, gs[1, 7:], summary)
    draw_field_panel(fig, gs[2, :])

    OUT.mkdir(parents=True, exist_ok=True)
    stem = OUT/"fig4_experimental_observability"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(stem)
    print("All symbols and sampled fields are theory placeholders; no experiment.")


if __name__ == "__main__":
    main()
