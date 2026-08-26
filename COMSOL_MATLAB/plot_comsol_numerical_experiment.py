"""Nature-style COMSOL scattering and finite-sample numerical experiment."""
from pathlib import Path
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "COMSOL_MATLAB" / "rounded_180k_results"
PER = DATA / "hard_boundary_scattering_scan"
FIN = DATA / "finite_20period_experiment"
OUT = ROOT / "figures" / "comsol_numerical_experiment"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 7.5,
    "axes.linewidth": 0.65,
    "xtick.major.width": 0.55,
    "ytick.major.width": 0.55,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "lines.linewidth": 1.15,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

BLUE = "#1764A2"
ORANGE = "#E67E22"
DARK = "#20262E"
GREY = "#77818C"


def read_csv(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}


coarse = read_csv(PER / "comsol_coarse_angle_frequency_map.csv")
fine = read_csv(PER / "comsol_fine_frequency_cuts.csv")
scanline = read_csv(FIN / "finite_bic_state_2_scanline.csv")
x = np.loadtxt(FIN / "finite_bic_field_x_m.csv", delimiter=",") * 1e3
y = np.loadtxt(FIN / "finite_bic_field_y_m.csv", delimiter=",") * 1e3
pr = np.loadtxt(FIN / "finite_bic_state_2_real_pa.csv", delimiter=",")
pi = np.loadtxt(FIN / "finite_bic_state_2_imag_pa.csv", delimiter=",")
pfield = pr + 1j * pi

fig = plt.figure(figsize=(7.2047244, 5.0393701))
gs = fig.add_gridspec(2, 12, height_ratios=[0.94, 1.06], hspace=0.43, wspace=1.65)
ax_a = fig.add_subplot(gs[0, 0:4])
ax_b = fig.add_subplot(gs[0, 4:8])
ax_c = fig.add_subplot(gs[0, 8:12])
ax_d = fig.add_subplot(gs[1, 0:8])
ax_e = fig.add_subplot(gs[1, 8:12])

# a, finite numerical experiment geometry.
ax_a.set_xlim(0, 148.4)
ax_a.set_ylim(-11, 67)
ax_a.add_patch(Rectangle((0, -10), 148.4, 10, fc="#59636E", ec="none"))
a = 7.42
for cell in range(20):
    x0 = cell * a
    ax_a.add_patch(Rectangle((x0 + 0.505, -5.15), 4.41, 5.15,
                             fc="white", ec="none"))
    ax_a.add_patch(Rectangle((x0 + 5.915, -1.32), 1.00, 1.32,
                             fc="white", ec="none"))
ax_a.plot([34.2, 114.2], [58, 58], color=DARK, lw=2.2, solid_capstyle="butt")
for xx in np.linspace(40, 108, 6):
    ax_a.annotate("", xy=(xx + 7.3, 40), xytext=(xx, 55),
                  arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=0.8))
ax_a.plot([0, 148.4], [20, 20], color=ORANGE, lw=0.8, ls=(0, (3, 2)))
ax_a.text(74.2, 61.2, "80-mm aperture", ha="center", va="bottom", color=DARK)
ax_a.text(146, 22.2, "scan line", ha="right", va="bottom", color=ORANGE)
ax_a.annotate("", xy=(148.4, -8.0), xytext=(0, -8.0),
              arrowprops=dict(arrowstyle="<->", color="white", lw=0.8))
ax_a.text(74.2, -7.4, "20 periods", ha="center", va="bottom", color="white")
ax_a.set_axis_off()

# b, global periodic COMSOL scattering landscape.
theta = np.unique(coarse["theta_deg"])
freq = np.unique(coarse["frequency_hz"]) / 1e3
assert np.all(np.diff(theta) > 0) and np.all(np.diff(freq) > 0)
z = np.empty((freq.size, theta.size))
for j, th in enumerate(theta):
    use = np.isclose(coarse["theta_deg"], th)
    order = np.argsort(coarse["frequency_hz"][use])
    z[:, j] = np.log10(np.maximum(coarse["large_groove_probe_abs_pa"][use][order], 1e-30))
mesh = ax_b.pcolormesh(theta, freq, z, shading="auto", cmap="magma")
ra = 1500 / (7.42e-3 * (1 + np.sin(np.deg2rad(theta)))) / 1e3
ax_b.plot(theta, ra, color="white", lw=1.0, ls=(0, (3, 2)))
ax_b.plot(7.0537, 180.046681327, "o", ms=3.8, mec="white", mew=0.7, mfc=ORANGE)
ax_b.set(xlabel=r"Bloch angle, $\theta$ (deg)", ylabel="Frequency (kHz)")
ax_b.set_xlim(theta.min(), theta.max())
ax_b.set_ylim(freq.min(), freq.max())
cb = fig.colorbar(mesh, ax=ax_b, fraction=0.047, pad=0.035)
cb.ax.set_title(r"$\log_{10}|p_{\rm cav}|$", fontsize=7.5, pad=3, loc="left")

# c, sub-hertz periodic-port signature at 7.10 degrees.
use = np.isclose(fine["theta_deg"], 7.10)
ff = fine["frequency_hz"][use]
order = np.argsort(ff)
ff = ff[order]
window = np.abs(ff - 180046.681327) <= 8
ff = ff[window]
eta0 = fine["eta_0"][use][order][window]
etam1 = fine["eta_m1"][use][order][window]
pc = fine["probe_abs_pa"][use][order][window]
pc = pc / np.nanmax(pc)
ax_c.plot(ff - 180046.681327, eta0, color=BLUE, label=r"$\eta_0$")
ax_c.plot(ff - 180046.681327, etam1, color=ORANGE, label=r"$\eta_{-1}$")
ax_c.plot(ff - 180046.681327, pc, color=DARK, ls=(0, (2, 1.5)), label=r"$|p_{\rm cav}|$")
ax_c.axvline(0, color="#BAC1C8", lw=0.65, zorder=0)
ax_c.set(xlabel=r"Detuning, $f-f_0$ (Hz)", ylim=(-0.03, 1.03))
ax_c.text(0.05, 0.86, r"$\theta=7.10^\circ$", transform=ax_c.transAxes,
          ha="left", va="top", color=GREY)
ax_c.legend(frameon=False, loc="upper right", handlelength=1.8, borderpad=0.1)

# d, finite-sample total field at the quasi-BIC frequency.
extent = [x.min(), x.max(), y.min(), y.max()]
real_field = np.real(pfield)
finite_values = np.abs(real_field[np.isfinite(real_field)])
lim = np.percentile(finite_values, 99.2)
im = ax_d.imshow(real_field, origin="lower", extent=extent, aspect="auto",
                 cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="bilinear")
ax_d.axhline(0, color=DARK, lw=0.7)
ax_d.axhline(20, color=ORANGE, lw=0.75, ls=(0, (3, 2)))
ax_d.set_xlim(0, 148.4)
ax_d.set_ylim(-6, 46)
ax_d.set(xlabel="Position, x (mm)", ylabel="Height, y (mm)")
ax_d.text(2.5, 41.5, r"$f=180.047$ kHz, $\theta=7.05^\circ$", color=DARK,
          ha="left", va="top")
cb2 = fig.colorbar(im, ax=ax_d, fraction=0.025, pad=0.018)
cb2.ax.set_title(r"Re$(p)$", fontsize=7.5, pad=3)

# e, directly measurable hydrophone trace at the dashed scan line in d.
xs = scanline["x_m"] * 1e3
ps = scanline["pressure_abs_pa"]
ax_e.axvspan(34.2, 114.2, color="#E9EDF1", lw=0, zorder=0)
ax_e.plot(xs, ps, color=ORANGE)
ax_e.set(xlabel="Position, x (mm)")
ax_e.set_title(r"$|p|$ at $y=20$ mm (Pa)", loc="left", pad=4)
ax_e.set_xlim(0, 148.4)
ax_e.text(74.2, 0.94, "source footprint", transform=ax_e.get_xaxis_transform(),
          ha="center", va="top", color=GREY)

panel_positions = [(-0.13, 1.06, DARK), (0.025, 0.975, "white"),
                   (0.025, 0.975, DARK), (-0.13, 1.06, DARK),
                   (0.025, 0.965, DARK)]
for label, ax, (xpos, ypos, color) in zip("abcde", [ax_a, ax_b, ax_c, ax_d, ax_e],
                                          panel_positions):
    ax.text(xpos, ypos, label, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="left", color=color, zorder=10)
for ax in [ax_b, ax_c, ax_d, ax_e]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.savefig(OUT / "Fig4_COMSOL_numerical_experiment.svg", bbox_inches="tight")
fig.savefig(OUT / "Fig4_COMSOL_numerical_experiment.pdf", bbox_inches="tight")
fig.savefig(OUT / "Fig4_COMSOL_numerical_experiment.png", dpi=600, bbox_inches="tight")
fig.savefig(OUT / "Fig4_COMSOL_numerical_experiment.tiff", dpi=600, bbox_inches="tight")
plt.close(fig)
