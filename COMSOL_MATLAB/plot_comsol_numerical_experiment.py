"""Nature-style numerical experiment for the finite Rayleigh-BIC sample."""
from pathlib import Path
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

ROOT = Path(__file__).resolve().parents[1]
POLE = ROOT / "Ni2019_MATLAB" / "results" / "fig4_wideangle_180k"
DRIVEN = ROOT / "Ni2019_MATLAB" / "results" / "fig4_experimental_observables_180k"
FIN = ROOT / "COMSOL_MATLAB" / "rounded_180k_results" / "finite_20period_experiment"
OUT = ROOT / "figures" / "comsol_numerical_experiment"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.labelsize": 8,
    "axes.titlesize": 8, "xtick.labelsize": 7.2, "ytick.labelsize": 7.2,
    "legend.fontsize": 7.2, "axes.linewidth": 0.65,
    "xtick.major.width": 0.55, "ytick.major.width": 0.55,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "lines.linewidth": 1.15, "svg.fonttype": "none", "pdf.fonttype": 42,
})

BLUE = "#1764A2"
ORANGE = "#E67E22"
DARK = "#20262E"
GREY = "#77818C"
LIGHT = "#DCE2E7"


def read_csv(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}


def angular_spectrum(path, frequency_hz):
    data = read_csv(path)
    x = data["x_m"]
    field = data["reflected_real_pa"] + 1j * data["reflected_imag_pa"]
    k0 = 2 * np.pi * frequency_hz / 1500
    # Evaluate the finite-window Fourier integral on a dense k grid. This is
    # a smooth display of the measured aperture spectrum; it does not improve
    # the physical resolution Delta kx ~ 2*pi/L.
    kx_over_k0 = np.linspace(-1.25, 1.25, 1201)
    phase = np.exp(-1j * np.outer(kx_over_k0 * k0, x))
    spectrum = phase @ (field * np.hanning(field.size))
    power = np.abs(spectrum) ** 2
    power_db = 10 * np.log10(np.maximum(power / np.max(power), 1e-6))
    return kx_over_k0, power_db


def measured_fwhm_from_complex_response(path):
    data = read_csv(path)
    theta_values = np.unique(data["theta_deg"])
    theta_out, width_out = [], []
    for angle in theta_values:
        sel = np.isclose(data["theta_deg"], angle, rtol=0, atol=2e-10)
        freq = data["frequency_hz"][sel]
        amp = data["Am1_real"][sel] + 1j * data["Am1_imag"][sel]
        valid = np.isfinite(freq) & np.isfinite(amp.real) & np.isfinite(amp.imag)
        if np.count_nonzero(valid) < 15:
            continue
        freq, amp = freq[valid], amp[valid]
        order = np.argsort(freq)
        freq, amp = freq[order], amp[order]
        intensity = np.abs(amp - 0.5 * (amp[0] + amp[-1])) ** 2
        peak_index = int(np.argmax(intensity))
        half = 0.5 * intensity[peak_index]
        left_candidates = np.flatnonzero(intensity[:peak_index] <= half)
        right_candidates = np.flatnonzero(intensity[peak_index + 1:] <= half)
        if left_candidates.size == 0 or right_candidates.size == 0:
            continue
        left_low = left_candidates[-1]
        right_high = peak_index + 1 + right_candidates[0]
        left = np.interp(half, intensity[left_low:left_low + 2],
                         freq[left_low:left_low + 2])
        # Reverse the descending segment so numpy.interp receives increasing xp.
        right = np.interp(half, intensity[right_high - 1:right_high + 1][::-1],
                          freq[right_high - 1:right_high + 1][::-1])
        if right > left:
            theta_out.append(angle)
            width_out.append(right - left)
    return np.array(theta_out), np.array(width_out)


pole = read_csv(POLE / "wideangle_physical_pole.csv")
theta = pole["theta_deg"]
fp = pole["frequency_hz"] / 1e3
linewidth = pole["linewidth_hz"]
assert np.all(np.diff(theta) > 0)

x = np.loadtxt(FIN / "finite_bic_field_x_m.csv", delimiter=",") * 1e3
y = np.loadtxt(FIN / "finite_bic_field_y_m.csv", delimiter=",") * 1e3
state_stems = [
    "gaussian_state_1_gamma",
    "gaussian_state_2_rayleigh_bic",
    "gaussian_state_3_above_bic",
]
state_angles = np.array([0.0, 7.05, 10.05404])
state_frequencies = np.array([178780.2005, 180046.681327, 180498.322])
state_titles = [r"$0^\circ$", r"$7.05^\circ$  (Rayleigh BIC)", r"$10.05^\circ$"]
incident_fields, reflected_fields = [], []
for stem in state_stems:
    incident_fields.append(
        np.loadtxt(FIN / f"{stem}_incident_real_pa.csv", delimiter=",")
        + 1j * np.loadtxt(FIN / f"{stem}_incident_imag_pa.csv", delimiter=",")
    )
    reflected_fields.append(
        np.loadtxt(FIN / f"{stem}_reflected_real_pa.csv", delimiter=",")
        + 1j * np.loadtxt(FIN / f"{stem}_reflected_imag_pa.csv", delimiter=",")
    )

theta_fwhm, driven_fwhm = measured_fwhm_from_complex_response(
    DRIVEN / "adaptive_scattering_response.csv")

# Figure contract: the finite scattering experiment exposes the opening of a
# second Floquet channel, while the projected resonant linewidth collapses at
# the same Rayleigh point as the outgoing eigenpole.
fig = plt.figure(figsize=(7.2047244, 7.05))
outer = fig.add_gridspec(3, 1, height_ratios=[0.86, 1.42, 0.78], hspace=0.42)
top = outer[0].subgridspec(1, 12, wspace=1.75)
ax_a = fig.add_subplot(top[0, 0:3])
ax_b = fig.add_subplot(top[0, 3:8])
ax_c = fig.add_subplot(top[0, 8:12])
field_grid = outer[1].subgridspec(2, 3, hspace=0.08, wspace=0.08)
field_axes = np.array([[fig.add_subplot(field_grid[r, c]) for c in range(3)]
                       for r in range(2)])
bottom = outer[2].subgridspec(1, 2, width_ratios=[1.25, 0.9], wspace=0.38)
ax_e1 = fig.add_subplot(bottom[0, 0])
ax_e2 = fig.add_subplot(bottom[0, 1])

# a, finite numerical experiment geometry.
ax_a.set_xlim(0, 148.4)
ax_a.set_ylim(-21, 67)
ax_a.add_patch(Rectangle((0, -10), 148.4, 10, fc="#59636E", ec="none"))
a = 7.42
for cell in range(20):
    x0 = cell * a
    ax_a.add_patch(Rectangle((x0 + 0.505, -5.15), 4.41, 5.15,
                             fc="white", ec="none"))
    ax_a.add_patch(Rectangle((x0 + 5.915, -1.32), 1.00, 1.32,
                             fc="white", ec="none"))
ax_a.plot([34.2, 114.2], [58, 58], color=DARK, lw=2.2, solid_capstyle="butt")
for xx in np.linspace(39, 107, 6):
    ax_a.annotate("", xy=(xx + 7.3, 39), xytext=(xx, 55),
                  arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=0.8))
ax_a.plot([0, 148.4], [20, 20], color=ORANGE, lw=0.8, ls=(0, (3, 2)))
ax_a.text(74.2, 61.2, "Gaussian aperture", ha="center", va="bottom", color=DARK)
ax_a.text(74.2, 22.2, "scan line", ha="center", va="bottom", color=ORANGE)
ax_a.annotate("", xy=(148.4, -14), xytext=(0, -14),
              arrowprops=dict(arrowstyle="<->", color=GREY, lw=0.8))
ax_a.text(74.2, -16.0, "20 periods", ha="center", va="top", color=GREY)
ax_a.set_axis_off()

# b, pole-derived resonance map: its bright linewidth pinches off at the BIC.
fgrid = np.linspace(178.55, 181.65, 1551)
spectral = np.zeros((fgrid.size, theta.size))
for j, (f0, gamma_hz) in enumerate(zip(fp, linewidth)):
    if gamma_hz <= 1e-12:
        continue
    half_khz = 0.5 * gamma_hz / 1e3
    spectral[:, j] = half_khz**2 / ((fgrid - f0)**2 + half_khz**2)
log_spectral = np.log10(np.maximum(spectral, 1e-6))
mesh = ax_b.pcolormesh(theta, fgrid, log_spectral, shading="nearest",
                       cmap="magma", vmin=-6, vmax=0, rasterized=True)
rayleigh = 1500 / (7.42e-3 * (1 + np.sin(np.deg2rad(theta)))) / 1e3
visible = (rayleigh >= fgrid.min()) & (rayleigh <= fgrid.max())
ax_b.plot(theta[visible], rayleigh[visible], color="white", lw=0.9,
          ls=(0, (3, 2)))
theta_bic = 7.0996629
ax_b.plot(theta_bic, 180.0, "o", ms=4.1, mec="white", mew=0.8, mfc=ORANGE)
ax_b.annotate("Rayleigh BIC", xy=(theta_bic, 180.0), xytext=(9.0, 179.35),
              color="white", fontsize=7,
              arrowprops=dict(arrowstyle="-", color="white", lw=0.7))
ax_b.set(xlabel=r"Bloch angle, $\theta$ (deg)", ylabel=r"$f$ (kHz)",
         xlim=(0, 20), ylim=(fgrid.min(), fgrid.max()))
ax_b.set_xticks([0, 5, 10, 15, 20])
cax = inset_axes(ax_b, width="28%", height="4%", loc="upper right", borderpad=1.0)
cb = fig.colorbar(mesh, cax=cax, orientation="horizontal", ticks=[-6, 0])
cb.ax.tick_params(labelsize=7.2, length=1.8, pad=1)
cb.ax.set_title(r"$\log_{10}\mathcal{L}$", fontsize=7.5, pad=2, loc="left")

# c, quantitative eigenpole linewidth collapse.
line_plot = np.maximum(linewidth, 1e-8)
ax_c.semilogy(theta, line_plot, color=ORANGE)
ax_c.axvline(theta_bic, color=LIGHT, lw=0.7, zorder=0)
ax_c.plot(theta_bic, 1e-8, "o", ms=4.2, mec=ORANGE, mew=0.9, mfc="white")
ax_c.set(xlabel=r"Bloch angle, $\theta$ (deg)", ylabel="Linewidth (Hz)",
         xlim=(-0.25, 20), ylim=(1e-8, 3e2))
ax_c.set_xticks([0, 5, 10, 15, 20])
ax_c.text(19.4, 130, "147 Hz", ha="right", va="bottom", color=GREY, fontsize=7)

# d, one shared physical scale for all six Gaussian-beam field maps.
use_y = (y >= 0) & (y <= 50)
all_values = []
for field in incident_fields + reflected_fields:
    values = np.real(field[use_y, :])
    all_values.append(np.abs(values[np.isfinite(values)]))
field_limit = np.percentile(np.concatenate(all_values), 99.4)
extent = [x.min(), x.max(), y[use_y].min(), y[use_y].max()]
field_image = None
for col in range(3):
    for row, fields in enumerate([incident_fields, reflected_fields]):
        ax = field_axes[row, col]
        field_image = ax.imshow(np.real(fields[col][use_y, :]), origin="lower",
                                extent=extent, aspect="auto", cmap="RdBu_r",
                                vmin=-field_limit, vmax=field_limit,
                                interpolation="bilinear", rasterized=True)
        ax.axhline(20, color=ORANGE, lw=0.65, ls=(0, (3, 2)))
        ax.axhline(0, color=DARK, lw=0.6)
        ax.set(xlim=(0, 148.4), ylim=(0, 50))
        ax.set_xticks([0, 50, 100, 150])
        ax.set_yticks([0, 20, 40])
        if row == 0:
            ax.set_title(state_titles[col], pad=3, color=DARK)
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("x (mm)")
        if col > 0:
            ax.tick_params(labelleft=False)
        else:
            ax.set_ylabel("y (mm)")
field_axes[0, 0].text(-0.32, 0.5, "Incident", rotation=90,
                      rotation_mode="anchor", transform=field_axes[0, 0].transAxes,
                      ha="center", va="center")
field_axes[1, 0].text(-0.32, 0.5, "Reflected", rotation=90,
                      rotation_mode="anchor", transform=field_axes[1, 0].transAxes,
                      ha="center", va="center")
field_cb = fig.colorbar(field_image, ax=field_axes.ravel().tolist(),
                        fraction=0.014, pad=0.012)
field_cb.ax.set_title(r"Re$(p)$ (Pa)", fontsize=7.2, pad=3)

# e(i), spatial Fourier content extracted from the complex reflected fields.
colors = [GREY, BLUE, ORANGE]
for stem, freq, color, label in zip(
        state_stems, state_frequencies, colors,
        [r"$0^\circ$", "Rayleigh BIC", r"$10.05^\circ$"]):
    kx_over_k0, power_db = angular_spectrum(FIN / f"{stem}_scanline.csv", freq)
    show = (kx_over_k0 >= -1.25) & (kx_over_k0 <= 1.25)
    ax_e1.plot(kx_over_k0[show], power_db[show], color=color, label=label)
ax_e1.axvspan(-1, 1, color="#F2F4F6", zorder=-3)
ax_e1.axvline(-1, color=LIGHT, lw=0.65)
ax_e1.axvline(1, color=LIGHT, lw=0.65)
ax_e1.text(0, -57, "propagating", ha="center", va="bottom", color=GREY, fontsize=7)
ax_e1.set(xlabel=r"Tangential spectrum, $k_x/k_0$", ylabel="Power (dB)",
          xlim=(-1.25, 1.25), ylim=(-60, 2))
ax_e1.set_xticks([-1, -0.5, 0, 0.5, 1])
theta_open = np.deg2rad(state_angles[-1])
g_over_k0 = 1500 / (7.42e-3 * state_frequencies[-1])
n0_position = -np.sin(theta_open)
nm1_position = g_over_k0 - np.sin(theta_open)
ax_e1.annotate(r"$n=0$", xy=(n0_position, -1.5), xytext=(-0.42, -18),
                fontsize=7, color=DARK,
                arrowprops=dict(arrowstyle="-", color=DARK, lw=0.65))
ax_e1.annotate(r"$n=-1$", xy=(nm1_position, -23), xytext=(0.68, -36),
                fontsize=7, color=ORANGE,
                arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.65))
ax_e1.legend(loc="upper right", ncol=3, handlelength=1.4,
             columnspacing=0.9, borderaxespad=0.2)

# e(ii), linewidth extracted from reflected fields overlays the eigenpole.
open_side = (theta >= theta_bic) & (theta <= 8.55)
ax_e2.semilogy(theta[open_side], np.maximum(linewidth[open_side], 1e-6),
               color=DARK, label="Eigenpole")
keep = (theta_fwhm >= theta_bic + 0.015) & (theta_fwhm <= 8.52)
indices = np.flatnonzero(keep)[::9]
ax_e2.semilogy(theta_fwhm[indices], driven_fwhm[indices], "o", ms=3.7,
               mfc="white", mec=ORANGE, mew=0.9, label="Driven scattering")
ax_e2.axvline(theta_bic, color=LIGHT, lw=0.7, zorder=-1)
ax_e2.set(xlabel=r"Angle, $\theta$ (deg)", ylabel="FWHM (Hz)",
          xlim=(7.08, 8.55), ylim=(1e-4, 8))
ax_e2.legend(loc="lower right", handlelength=1.5)

# Panel labels and common styling.
ax_a.text(-0.13, 1.06, "a", transform=ax_a.transAxes, fontsize=9,
          fontweight="bold", va="top", color=DARK)
ax_b.text(0.02, 0.97, "b", transform=ax_b.transAxes, fontsize=9,
          fontweight="bold", va="top", color="white", zorder=10)
ax_c.text(0.03, 0.97, "c", transform=ax_c.transAxes, fontsize=9,
          fontweight="bold", va="top", color=DARK)
field_axes[0, 0].text(-0.20, 1.16, "d", transform=field_axes[0, 0].transAxes,
                      fontsize=9, fontweight="bold", va="top", color=DARK)
ax_e1.text(-0.13, 1.06, "e", transform=ax_e1.transAxes, fontsize=9,
           fontweight="bold", va="top", color=DARK)

for ax in [ax_b, ax_c, ax_e1, ax_e2, *field_axes.ravel()]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.savefig(OUT / "Fig4_COMSOL_numerical_experiment.svg", bbox_inches="tight")
fig.savefig(OUT / "Fig4_COMSOL_numerical_experiment.pdf", bbox_inches="tight")
fig.savefig(OUT / "Fig4_COMSOL_numerical_experiment.png", dpi=600, bbox_inches="tight")
fig.savefig(OUT / "Fig4_COMSOL_numerical_experiment.tiff", dpi=600, bbox_inches="tight")
plt.close(fig)
