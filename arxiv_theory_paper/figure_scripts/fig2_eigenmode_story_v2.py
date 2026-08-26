"""Fig. 2 | Eigenvalue and radiation certification of the 180-kHz Rayleigh BIC.

All numerical panels use the same 180-kHz, 7.10-degree geometry and the same
continued complex-eigenfrequency branch.  The evidence chain is deliberately
direct: Re(f_p) crosses the Rayleigh line, Im(f_p) vanishes, and both exterior
radiation amplitudes vanish at the same Bloch angle.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from scipy.io import loadmat


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7.8,
    "axes.labelsize": 8.0,
    "axes.linewidth": 0.75,
    "xtick.labelsize": 7.4,
    "ytick.labelsize": 7.4,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 2.7,
    "ytick.major.size": 2.7,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "mathtext.fontset": "dejavusans",
})

INK = "#20252B"
NAVY = "#164F86"
BLUE = "#2675A8"
ORANGE = "#DC7A16"
RED = "#D33E32"
GRAY = "#7D8790"
LIGHT = "#D9DEE2"
SLAB = "#B9C8CF"
SLAB_DARK = "#758B98"
MODE = "#A7D8E8"
CMAP = LinearSegmentedColormap.from_list(
    "rayleigh_pole", ["#101C35", "#164B68", "#168382", "#79BF7A", "#F3D36C"]
)

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "arxiv_theory_paper" / "figures"
GEOM_MAT = ROOT / "Ni2019_MATLAB" / "results" / "StrictRayleighBIC_180kHz_7p10deg_final.mat"
EVOL_CSV = ROOT / "Ni2019_MATLAB" / "results" / "kx_channel_evolution_180k" / "kx_channel_evolution.csv"
FIELD_MAT = ROOT / "Ni2019_MATLAB" / "results" / "fig2_eigenmode_fields_180k" / "eigenmode_fields.mat"
C0 = 1500.0


def panel_label(ax, letter, inside=False):
    x, y = ((0.018, 0.972) if inside else (-0.09, 1.06))
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9.2,
            fontweight="bold", ha="left", va="top",
            color="white" if inside else INK, clip_on=False, zorder=30)


def clean(ax, grid=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(color=LIGHT, lw=0.55, alpha=0.65)
        ax.set_axisbelow(True)


def load_geometry():
    d = loadmat(GEOM_MAT, squeeze_me=True, struct_as_record=False)
    x = np.asarray(d["x"], dtype=float)
    depths = np.asarray(d["depths"], dtype=float)
    dims = np.asarray(d["dims"], dtype=float)
    a_mm = float(d["a"]) * 1e3
    return {
        "kappa": x[0], "Omega": 1.0 - x[0],
        "theta": np.degrees(np.arcsin(x[0] / (1.0 - x[0]))),
        "a_mm": a_mm, "d1_mm": depths[0], "d2_mm": depths[1],
        "w1_mm": dims[0], "w2_mm": dims[1], "gap_mm": dims[2],
    }


def dim_arrow(ax, p0, p1, text, offset=(0, 0)):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="<->", mutation_scale=6,
                                 lw=0.7, color=GRAY))
    xm = (p0[0] + p1[0]) / 2 + offset[0]
    ym = (p0[1] + p1[1]) / 2 + offset[1]
    ax.text(xm, ym, text, ha="center", va="center", color=INK, fontsize=7.3)


def draw_structure(ax, g):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    top = np.array([[0.04, 0.71], [0.83, 0.71], [0.96, 0.82], [0.17, 0.82]])
    front = np.array([[0.04, 0.71], [0.83, 0.71], [0.83, 0.66], [0.04, 0.66]])
    side = np.array([[0.83, 0.71], [0.96, 0.82], [0.96, 0.77], [0.83, 0.66]])
    ax.add_patch(Polygon(front, facecolor=SLAB_DARK, edgecolor=INK, lw=0.65))
    ax.add_patch(Polygon(side, facecolor="#8796AA", edgecolor=INK, lw=0.65))
    ax.add_patch(Polygon(top, facecolor=SLAB, edgecolor=INK, lw=0.7))
    vx, vz = np.array([0.13, 0.0]), np.array([0.13, 0.11])
    for j in range(6):
        origin = np.array([0.055, 0.718]) + j * vx
        if j:
            q = origin - 0.010 * vz
            ax.plot([q[0], q[0]+vz[0]], [q[1], q[1]+vz[1]], color="white", lw=0.55)
        for frac, width in [(0.12, 0.064), (0.73, 0.026)]:
            p = origin + frac*vx + 0.22*vz
            aperture = np.array([p, p + width/0.13*vx,
                                 p + width/0.13*vx + 0.52*vz, p + 0.52*vz])
            ax.add_patch(Polygon(aperture, facecolor="white", edgecolor=INK, lw=0.45))
    ax.text(0.50, 0.91, "periodic acoustic metasurface", ha="center", fontsize=7.3, color=GRAY)

    x0, y0, width, height = 0.12, 0.19, 0.76, 0.25
    ax.add_patch(Rectangle((x0, y0), width, height, facecolor=SLAB, edgecolor=INK, lw=0.75))
    w1, w2 = g["w1_mm"]/g["a_mm"], g["w2_mm"]/g["a_mm"]
    gap = g["gap_mm"]/g["a_mm"]
    d1, d2 = g["d1_mm"]/g["a_mm"], g["d2_mm"]/g["a_mm"]
    margin = (1-w1-w2-gap)/2
    gx1, gx2 = x0 + width*margin, x0 + width*(margin+w1+gap)
    ax.add_patch(Rectangle((gx1, y0+height*(1-d1)), width*w1, height*d1,
                           facecolor="white", edgecolor=INK, lw=0.7))
    ax.add_patch(Rectangle((gx2, y0+height*(1-d2)), width*w2, height*d2,
                           facecolor="white", edgecolor=INK, lw=0.7))
    ax.plot([x0, x0+width], [y0+height, y0+height], color=INK, lw=0.85)
    dim_arrow(ax, (x0, 0.12), (x0+width, 0.12), "a", offset=(0, -0.027))
    ax.text(gx1+width*w1/2, y0+height+0.035, r"$w_1$", ha="center", fontsize=7.3)
    ax.text(gx2+width*w2/2, y0+height+0.035, r"$w_2$", ha="center", fontsize=7.3)
    ax.text(gx1+width*w1+width*gap/2, y0+height+0.035, "g", ha="center", fontsize=7.3)
    ax.text(gx1-0.018, y0+height*(1-d1/2), r"$d_1$", ha="right", va="center", fontsize=7.3)
    ax.text(gx2+width*w2+0.014, y0+height*(1-d2/2), r"$d_2$", ha="left", va="center", fontsize=7.3)


def load_evolution():
    d = np.genfromtxt(EVOL_CSV, delimiter=",", names=True)
    return d[np.argsort(d["kappa"])]


def branch_data(e, g):
    theta = np.degrees(np.arcsin(e["kappa"] / e["Omega_real"]))
    factor = C0 / (g["a_mm"]*1e-3) / 1e3
    fp = e["Omega_real"] * factor
    gamma = np.abs(e["Omega_imag"]) * factor
    return theta, fp, gamma


def draw_dispersion(ax, e, g):
    theta, fp, gamma = branch_data(e, g)
    span = min(g["theta"] - theta.min(), theta.max() - g["theta"])
    tg = g["theta"] + np.linspace(-span, span, 501)
    fp_i = np.interp(tg, theta, fp)
    gamma_i = np.interp(tg, theta, gamma)
    fg = 180.0 + np.linspace(-2.45, 2.45, 501)
    pole_distance = np.sqrt((fg[:, None] - fp_i[None, :])**2 + gamma_i[None, :]**2)
    pole_distance = np.maximum(pole_distance, 1e-10)
    proximity = -np.log10(pole_distance / 180.0)
    im = ax.imshow(proximity, origin="lower", aspect="auto",
                   extent=[tg.min(), tg.max(), fg.min(), fg.max()],
                   interpolation="nearest", cmap=CMAP, vmin=1.8, vmax=10.0)
    fra = C0 / ((g["a_mm"]*1e-3)*(1+np.sin(np.deg2rad(tg)))) / 1e3
    ax.plot(tg, fra, color=ORANGE, lw=1.25, ls=(0, (4, 2.6)))
    ax.plot(g["theta"], 180.0, "o", ms=6.0, mfc="none", mec=RED, mew=1.2, zorder=8)
    ax.annotate(r"$\mathrm{Re}\,f_p=f_{\rm RA}$", (g["theta"], 180.0),
                xytext=(g["theta"]+0.34, 180.68), color="white", fontsize=7.3,
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.7))
    ax.text(g["theta"]+0.91, np.interp(g["theta"]+0.91, tg, fp_i)+0.15,
            r"$\mathrm{Re}\,f_p$", color="white", fontsize=7.3, ha="center")
    ax.text(6.35, np.interp(6.35, tg, fra)+0.18, r"$f_{\rm RA}$",
            color=ORANGE, fontsize=7.3, rotation=-39, rotation_mode="anchor")
    ax.set_xlabel("Bloch angle θ (deg)")
    ax.set_ylabel("Frequency (kHz)")
    ax.set_xlim(tg.min(), tg.max())
    ax.set_ylim(fg.min(), fg.max())
    clean(ax)
    cax = inset_axes(ax, width="2.4%", height="82%", loc="lower left",
                     bbox_to_anchor=(1.018, 0.09, 1, 1), bbox_transform=ax.transAxes,
                     borderpad=0)
    cb = plt.colorbar(im, cax=cax, ticks=[2, 6, 10])
    cb.set_label(r"$-\log_{10}(|f-f_p|/f_{\rm BIC})$", labelpad=2)
    cb.ax.tick_params(length=2.0, labelsize=7.3)


def draw_complex_eigenvalue(ax, e, g):
    """Real and imaginary parts of one continued eigenvalue branch."""
    theta, fp, gamma_khz = branch_data(e, g)
    gamma_hz = 1e3 * gamma_khz
    tg = np.linspace(theta.min(), theta.max(), 501)
    fra = C0 / ((g["a_mm"]*1e-3)*(1+np.sin(np.deg2rad(tg)))) / 1e3

    line_re, = ax.plot(theta, fp, color=NAVY, lw=1.55, label=r"$\mathrm{Re}\,f_p$")
    line_ra, = ax.plot(tg, fra, color=ORANGE, lw=1.2, ls=(0, (4, 2.5)),
                       label=r"$f_{\rm RA}$")
    ax.axvline(g["theta"], color=GRAY, lw=0.75, ls=(0, (2.2, 2.2)), zorder=0)
    ax.plot(g["theta"], 180.0, "o", ms=5.6, mfc=RED, mec="white", mew=0.75, zorder=8)
    ax.annotate(r"$\mathrm{Re}\,f_p=f_{\rm RA}$", (g["theta"], 180.0),
                xytext=(g["theta"]+0.32, 180.56), color=RED, fontsize=7.3,
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.7))
    ax.set_xlabel("Bloch angle θ (deg)")
    ax.set_ylabel("Frequency (kHz)")
    ax.set_xlim(theta.min(), theta.max())
    ax.set_ylim(177.4, 182.5)
    clean(ax, grid=True)

    axr = ax.twinx()
    root = np.argmin(np.abs(e["delta_kappa"]))
    for sl in (slice(0, root), slice(root+1, None)):
        axr.semilogy(theta[sl], gamma_hz[sl], "o-", color=RED, lw=1.0,
                     ms=2.6, mfc="white", mec=RED, mew=0.65)
    floor = 1.1e-8
    axr.plot(g["theta"], floor, "o", ms=5.0, mfc=RED, mec="white", mew=0.75,
             clip_on=False, zorder=8)
    axr.annotate(r"$|\mathrm{Im}\,f_p|=0$", xy=(g["theta"], floor),
                 xytext=(g["theta"]+0.43, 5e-8), color=RED, fontsize=7.3,
                 arrowprops=dict(arrowstyle="->", color=RED, lw=0.7))
    axr.set_ylabel(r"$|\mathrm{Im}\,f_p|$ (Hz)", color=RED)
    axr.set_ylim(8e-9, 5)
    axr.tick_params(axis="y", colors=RED)
    axr.spines["right"].set_color(RED)
    axr.spines["top"].set_visible(False)
    ax.legend(handles=[line_re,line_ra], loc="lower left", frameon=False,
              fontsize=7.3, handlelength=1.55, borderaxespad=0.4)


def draw_imaginary_part(ax, e, g):
    """Radiative decay rate from the imaginary part of the same eigenvalue."""
    theta, _, gamma_khz = branch_data(e, g)
    gamma_hz = 1e3 * gamma_khz
    dk = e["delta_kappa"]
    for side, color, fill, label in [(-1, NAVY, NAVY, r"$\theta<\theta_{\rm BIC}$"),
                                     (1, ORANGE, "white", r"$\theta>\theta_{\rm BIC}$")]:
        m = (np.sign(dk) == side) & (gamma_hz > 0)
        ax.semilogy(theta[m], gamma_hz[m], "o-", color=color, lw=1.15,
                    ms=2.8, mfc=fill, mec=color, mew=0.65, label=label)
    floor = 1.2e-8
    ax.axvline(g["theta"], color=RED, lw=0.75, ls=(0, (2.2, 2.2)))
    ax.plot(g["theta"], floor, "o", ms=5.2, mfc=RED, mec="white", mew=0.75,
            clip_on=False, zorder=8)
    ax.annotate(r"$\mathrm{Im}\,f_p=0$", xy=(g["theta"], floor),
                xytext=(g["theta"]+0.34, 8e-8), color=RED, fontsize=7.3,
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.7))
    ax.set_xlabel("Bloch angle θ (deg)")
    ax.set_ylabel(r"$|\mathrm{Im}\,f_p|$ (Hz)")
    ax.set_xlim(theta.min(), theta.max())
    ax.set_ylim(8e-9, 5)
    ax.legend(loc="upper left", fontsize=7.3, frameon=False,
              handlelength=1.45, borderaxespad=0.35)
    clean(ax, grid=True)


def draw_linewidth_cloud(ax, e, g):
    """Single-pole radiative contrast; its width and contrast vanish at BIC."""
    theta, _, gamma = branch_data(e, g)
    span = 0.52
    tg = g["theta"] + np.linspace(-span, span, 401)
    # ``branch_data`` returns kHz, whereas the local detuning axis is in Hz.
    gg = 1e3 * np.interp(tg, theta, gamma)
    detuning = np.linspace(-0.65, 0.65, 401)  # Hz relative to Re f_p(theta)
    numerator = gg[None, :] ** 2
    denominator = detuning[:, None] ** 2 + numerator
    contrast = np.divide(numerator, denominator, out=np.zeros_like(denominator),
                         where=denominator > 0)
    safe_contrast = np.maximum(contrast, 1e-7)
    assert np.all(safe_contrast > 0)
    log_contrast = np.log10(safe_contrast)
    im = ax.imshow(log_contrast, origin="lower", aspect="auto",
                   extent=[tg.min(), tg.max(), detuning.min(), detuning.max()],
                   interpolation="nearest", cmap=CMAP, vmin=-7, vmax=0)
    ax.axvline(g["theta"], color="white", lw=0.8, ls=(0, (2.2, 2.2)))
    ax.plot(g["theta"], 0, "o", ms=5.0, mfc=RED, mec="white", mew=0.75, zorder=8)
    ax.text(g["theta"]+0.055, 0.48, "vanishing linewidth", color="white",
            fontsize=7.3, ha="left", va="center")
    ax.set_xlabel("Bloch angle θ (deg)")
    ax.set_ylabel(r"$f-\mathrm{Re}\,f_p$ (Hz)")
    clean(ax)
    cax = inset_axes(ax, width="31%", height="3.8%", loc="lower left",
                     bbox_to_anchor=(0.055, 0.075, 1, 1), bbox_transform=ax.transAxes,
                     borderpad=0)
    cb = plt.colorbar(im, cax=cax, orientation="horizontal", ticks=[-6, -3, 0])
    cb.ax.tick_params(length=1.8, labelsize=7.0, colors="white", pad=1.2)
    cb.outline.set_edgecolor("white")
    cb.ax.text(1.12, 0.5, r"$\log_{10}\mathcal{R}_{\rm pole}$",
               transform=cb.ax.transAxes, color="white", fontsize=7.3,
               ha="left", va="center")


def draw_q(ax, e):
    dk = e["delta_kappa"]
    valid = (dk != 0) & (np.abs(e["Omega_imag"]) > 1e-15)
    q = np.full(e.shape, np.nan, dtype=float)
    q[valid] = np.abs(e["Omega_real"][valid]) / (2*np.abs(e["Omega_imag"][valid]))
    for side, color, fill, label in [(-1, NAVY, NAVY, "below opening"),
                                     (1, ORANGE, "white", "above opening")]:
        m = valid & (np.sign(dk) == side)
        x, y = np.abs(dk[m]), q[m]
        o = np.argsort(x)
        ax.loglog(x[o], y[o], "o-", ms=3.0, lw=1.05, color=color,
                  mfc=fill, mec=color, mew=0.7, label=label)
    ax.annotate("Q → ∞", xy=(1.1e-6, 1.1e12), xytext=(5.0e-6, 1.7e12),
                color=RED, fontsize=7.3,
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.7))
    ax.set_xlabel(r"$|\kappa-\kappa_{\mathrm{BIC}}|$")
    ax.set_ylabel("")
    ax.text(0.025, 0.51, "Q", transform=ax.transAxes, rotation=90,
            rotation_mode="anchor",
            ha="left", va="center", color=INK, fontsize=8.0)
    ax.set_xlim(8e-7, 3.2e-2)
    ax.set_ylim(1e4, 5e12)
    ax.legend(loc="lower left", fontsize=7.3, frameon=False, handlelength=1.5)
    clean(ax, grid=True)


def draw_channel_zero(ax, e, g):
    theta = np.degrees(np.arcsin(e["kappa"] / e["Omega_real"]))
    dk = e["delta_kappa"]
    for values, color, label in [(e["A0_over_C"], BLUE, r"$|A_0|/\Vert C\Vert_2$"),
                                  (e["Am1_over_C"], ORANGE, r"$|A_{-1}|/\Vert C\Vert_2$")]:
        for m in (dk < 0, dk > 0):
            ax.semilogy(theta[m], values[m], color=color, lw=1.25)
        ax.text(theta[-5], values[-5]*1.25, label, color=color, fontsize=7.3,
                ha="right", va="bottom")
    floor = 1.1e-7
    ax.plot(g["theta"], floor, "o", ms=5.2, mfc=RED, mec="white", mew=0.75,
            clip_on=False, zorder=8)
    ax.axvline(g["theta"], color=RED, lw=0.75, ls=(0, (2.2, 2.2)))
    ax.annotate(r"$A_0=A_{-1}=0$", xy=(g["theta"], floor),
                xytext=(g["theta"]+0.55, 3.4e-7), color=RED, fontsize=7.3,
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.7))
    ax.set_xlabel("Bloch angle θ (deg)")
    ax.set_ylabel("Normalized radiation amplitude")
    ax.set_xlim(theta.min(), theta.max())
    ax.set_ylim(8e-8, 2e-2)
    clean(ax, grid=True)


def load_eigenfields():
    d = loadmat(FIELD_MAT, squeeze_me=True, struct_as_record=False)
    return d["x"], d["y"], list(np.atleast_1d(d["states"]))


def draw_bic_eigenfield(ax, fig, g):
    x, y, states = load_eigenfields()
    state = states[1]
    # |p| is Bloch-periodic at the real-frequency BIC, so two repeated cells
    # show the periodic localization without introducing an incident field.
    z = np.asarray(state.log_magnitude)
    z2 = np.concatenate([z[:, :-1], z], axis=1)
    ax.set_facecolor("#D6D9DC")
    im = ax.imshow(z2, origin="lower", extent=[0, 2, y.min(), y.max()],
                   aspect="equal", interpolation="nearest", cmap=CMAP,
                   vmin=-5.5, vmax=0)
    ax.axhline(0, color=INK, lw=0.7)
    margin = (1-g["w1_mm"]/g["a_mm"]-g["w2_mm"]/g["a_mm"]-g["gap_mm"]/g["a_mm"])/2
    w1, w2 = g["w1_mm"]/g["a_mm"], g["w2_mm"]/g["a_mm"]
    gap = g["gap_mm"]/g["a_mm"]
    d1, d2 = g["d1_mm"]/g["a_mm"], g["d2_mm"]/g["a_mm"]
    for cell in (0, 1):
        ax.add_patch(Rectangle((cell+margin, -d1), w1, d1,
                               fill=False, edgecolor="white", lw=0.55))
        ax.add_patch(Rectangle((cell+margin+w1+gap, -d2), w2, d2,
                               fill=False, edgecolor="white", lw=0.55))
    ax.set_xlim(0, 2); ax.set_ylim(-0.74, 1.05)
    ax.set_xticks([0, 1, 2]); ax.set_yticks([-0.5, 0, 0.5, 1.0])
    ax.set_xlabel(r"$x/a$"); ax.set_ylabel(r"$y/a$")
    ax.set_title(r"Rayleigh BIC eigenfield, $\theta=7.10^\circ$",
                 fontsize=7.3, color=INK, pad=3.0)
    for spine in ax.spines.values():
        spine.set_linewidth(0.65)
    panel_label(ax, "d")
    cax = inset_axes(ax, width="3.6%", height="88%", loc="lower left",
                     bbox_to_anchor=(1.035, 0.06, 1, 1), bbox_transform=ax.transAxes,
                     borderpad=0)
    cb = fig.colorbar(im, cax=cax, ticks=[-5, -3, -1, 0])
    cb.set_label(r"$\log_{10}|p_{\rm eig}/p_{\max}|$", labelpad=2)
    cb.ax.tick_params(length=2.0, labelsize=7.3)


def main():
    g, e = load_geometry(), load_evolution()
    root = np.argmin(np.abs(e["delta_kappa"]))
    assert np.isclose(g["kappa"], e["kappa"][root])
    assert e["A0_over_C"][root] == 0 and e["Am1_over_C"][root] == 0
    fig = plt.figure(figsize=(7.2047, 4.78), facecolor="white")
    gs = fig.add_gridspec(2, 18, height_ratios=[1.04, 0.96],
                          left=0.055, right=0.955, bottom=0.095, top=0.965,
                          wspace=1.70, hspace=0.40)
    axa = fig.add_subplot(gs[0, :6])
    axb = fig.add_subplot(gs[0, 6:])
    axc = fig.add_subplot(gs[1, :8])
    axd = fig.add_subplot(gs[1, 10:])
    draw_structure(axa, g)
    draw_complex_eigenvalue(axb, e, g)
    draw_channel_zero(axc, e, g)
    draw_bic_eigenfield(axd, fig, g)
    panel_label(axa, "a")
    panel_label(axb, "b")
    panel_label(axc, "c")
    FIG.mkdir(parents=True, exist_ok=True)
    stem = FIG / "fig2_eigenmode_story_v2"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.025)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.025)
    fig.savefig(stem.with_suffix(".png"), dpi=450, bbox_inches="tight", pad_inches=0.025)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pad_inches=0.025)
    plt.close(fig)


if __name__ == "__main__":
    main()
