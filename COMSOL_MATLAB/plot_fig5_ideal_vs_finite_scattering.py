"""Figure 5: finite Gaussian beams regularize an ideal Rayleigh BIC."""
from pathlib import Path
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, Circle

ROOT = Path(__file__).resolve().parents[1]
POLE = ROOT / "Ni2019_MATLAB" / "results" / "fig4_wideangle_180k"
SCAT = (ROOT / "Ni2019_MATLAB" / "results"
        / "fig4_experimental_observables_180k")
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

DARK = "#20262E"
GREY = "#7A848E"
LIGHT = "#DCE3E8"
BLUE = "#1764A2"
CYAN = "#62AFC5"
ORANGE = "#E67E22"
PURPLE = "#7356A5"
THETA_BIC = 7.09966291282572
F_BIC = 180000.0
C0 = 1500.0
F0 = 180e3
PERIOD_MM = 7.42
CURRENT_APERTURE_MM = 80.0
GAUSSIAN_FWHM_FACTOR = 2.3548200450309493


def read_csv(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}


def gaussian_weight(theta_internal, theta_center, waist_mm):
    """Normalized angular intensity for exp[-(x/w)^2] pressure amplitude."""
    k0 = 2 * np.pi * F0 / C0
    delta_kx = k0 * (np.sin(np.deg2rad(theta_internal))
                     - np.sin(np.deg2rad(theta_center)))
    weight = np.exp(-0.5 * (waist_mm * 1e-3 * delta_kx) ** 2)
    return weight / np.trapezoid(weight, theta_internal)


def spectral_model(pole):
    """Pole spectral density and its finite-beam angular convolutions."""
    assert np.all(np.diff(pole["theta_deg"]) > 0)
    theta_internal = np.linspace(-20, 20, 2001)
    pole_frequency = np.interp(
        np.abs(theta_internal), pole["theta_deg"], pole["frequency_hz"])
    linewidth = np.interp(
        np.abs(theta_internal), pole["theta_deg"], pole["linewidth_hz"])
    frequency = np.arange(178600.0, 181402.0, 2.0)
    half = 0.5 * linewidth
    ideal_internal = half[None, :] ** 2 / (
        (frequency[:, None] - pole_frequency[None, :]) ** 2
        + half[None, :] ** 2 + 1e-300)

    theta_command = np.linspace(5.8, 8.5, 136)
    pole_command = np.interp(
        theta_command, pole["theta_deg"], pole["frequency_hz"])
    linewidth_command = np.interp(
        theta_command, pole["theta_deg"], pole["linewidth_hz"])
    half_command = 0.5 * linewidth_command
    ideal_command = half_command[None, :] ** 2 / (
        (frequency[:, None] - pole_command[None, :]) ** 2
        + half_command[None, :] ** 2 + 1e-300)

    waist_grid = np.array([25., 35., 50., 75., 100., 150., 200., 250., 300.])
    maps, angular_bandwidths = {}, {}
    dtheta = theta_internal[1] - theta_internal[0]
    trap = np.ones(theta_internal.size)
    trap[[0, -1]] = 0.5
    for waist in waist_grid:
        weights = np.column_stack([
            gaussian_weight(theta_internal, theta, waist)
            for theta in theta_command])
        response = ideal_internal @ (weights * trap[:, None] * dtheta)
        maps[waist] = response
        integration = weights * trap[:, None] * dtheta
        mean_frequency = np.sum(pole_frequency[:, None] * integration, axis=0)
        variance = np.sum(
            (pole_frequency[:, None] - mean_frequency[None, :]) ** 2
            * integration, axis=0)
        angular_bandwidths[waist] = GAUSSIAN_FWHM_FACTOR * np.sqrt(variance)
    return (theta_command, frequency, ideal_command, maps, angular_bandwidths,
            waist_grid, linewidth_command, pole_command)


pole = read_csv(POLE / "wideangle_physical_pole.csv")
radiation = read_csv(SCAT / "pole_and_radiation_summary.csv")
(theta, frequency, ideal_map, gaussian_maps, angular_bandwidths,
 waist_grid, eigen_width, pole_frequency) = spectral_model(pole)

# Source data exported from the exact curves drawn in quantitative panels.
with open(OUT / "Fig5_gaussian_beam_scaling.csv", "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["theta_deg", "intrinsic_linewidth_hz",
                     "angular_bandwidth_w25_hz", "angular_bandwidth_w100_hz",
                     "angular_bandwidth_w200_hz"])
    writer.writerows(zip(theta, eigen_width, angular_bandwidths[25.],
                         angular_bandwidths[100.], angular_bandwidths[200.]))
with open(OUT / "Fig5_minimum_width_vs_waist.csv", "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["waist_mm", "angular_fwhm_deg",
                     "minimum_angular_bandwidth_hz", "matched_aperture_mm",
                     "minimum_period_count"])
    k0 = 2 * np.pi * F0 / C0
    for waist in waist_grid:
        angular_fwhm = np.rad2deg(
            GAUSSIAN_FWHM_FACTOR / (k0 * waist * 1e-3
                                    * np.cos(np.deg2rad(THETA_BIC))))
        aperture = 3.2 * waist
        writer.writerow([waist, angular_fwhm,
                         np.nanmin(angular_bandwidths[waist]), aperture,
                         np.ceil(aperture / PERIOD_MM)])

# Figure contract: a finite beam cannot itself realize the single-k BIC limit;
# its known angular spectrum permits a controlled extrapolation to that limit.
fig = plt.figure(figsize=(7.2047244, 7.15))
outer = fig.add_gridspec(3, 1, height_ratios=[0.78, 1.06, 0.88], hspace=0.47)
top = outer[0].subgridspec(1, 2, width_ratios=[1.18, 0.82], wspace=0.28)
middle = outer[1].subgridspec(1, 2, wspace=0.16)
bottom = outer[2].subgridspec(1, 3, wspace=0.40)
ax_a = fig.add_subplot(top[0, 0])
ax_b = fig.add_subplot(top[0, 1])
ax_c = fig.add_subplot(middle[0, 0])
ax_d = fig.add_subplot(middle[0, 1])
ax_e = fig.add_subplot(bottom[0, 0])
ax_f = fig.add_subplot(bottom[0, 1])
ax_g = fig.add_subplot(bottom[0, 2])

# a, measurement and k-space recovery principle.
ax_a.set(xlim=(0, 1), ylim=(0, 1))
ax_a.axis("off")
sample_left, sample_right, sample_y = 0.04, 0.70, 0.19
ax_a.add_patch(Rectangle((sample_left, sample_y - 0.055),
                         sample_right - sample_left, 0.055,
                         facecolor="#59636D", edgecolor="none"))
cell_width = (sample_right - sample_left) / 20
for xg in np.linspace(sample_left, sample_right - cell_width, 20):
    ax_a.add_patch(Rectangle((xg + 0.18 * cell_width, sample_y - 0.035),
                             0.28 * cell_width, 0.035,
                             facecolor="white", edgecolor="none"))
    ax_a.add_patch(Rectangle((xg + 0.68 * cell_width, sample_y - 0.025),
                             0.13 * cell_width, 0.025,
                             facecolor="white", edgecolor="none"))
ax_a.text(0.37, 0.085, r"20 periods ($20a=148.4$ mm)", ha="center")
source_l, source_r, source_y = 0.20, 0.56, 0.87
ax_a.plot([source_l, source_r], [source_y, source_y], color=ORANGE,
          lw=4.0, solid_capstyle="butt")
ax_a.annotate("", xy=(source_l, 0.94), xytext=(source_r, 0.94),
              arrowprops=dict(arrowstyle="<->", color=DARK, lw=0.75))
ax_a.text(0.38, 0.965, "80-mm active aperture", ha="center", va="bottom")
beam = Polygon([[source_l, source_y - 0.015], [source_r, source_y - 0.015],
                [0.54, sample_y], [0.21, sample_y]], closed=True,
               facecolor="#E7F2F5", edgecolor="none", alpha=0.95)
ax_a.add_patch(beam)
for offset, alpha in [(-0.10, 0.25), (-0.05, 0.48), (0, 0.85),
                      (0.05, 0.48), (0.10, 0.25)]:
    ax_a.plot([0.38 + offset, 0.38 + offset - 0.04],
              [source_y - 0.02, sample_y + 0.02], color=CYAN,
              lw=0.9, alpha=alpha)
ax_a.text(0.39, 0.67, r"$w=25$ mm", ha="center", color=ORANGE)
scan_y = 0.43
ax_a.plot([sample_left, sample_right], [scan_y, scan_y], color=BLUE,
          lw=0.8, ls=(0, (3, 2)))
for xp in np.linspace(0.10, 0.64, 6):
    ax_a.add_patch(Circle((xp, scan_y), 0.010, facecolor="white",
                          edgecolor=BLUE, lw=0.7))
ax_a.text(0.37, scan_y + 0.035, r"complex scan $p_{\rm refl}(x,f)$",
          ha="center", color=BLUE)
ax_a.annotate("", xy=(0.82, 0.51), xytext=(0.70, 0.51),
              arrowprops=dict(arrowstyle="->", color=DARK, lw=0.8))
ax_a.text(0.76, 0.55, r"$\mathcal{F}_x$", ha="center")
ax_a.plot([0.84, 0.96], [0.51, 0.51], color=DARK, lw=0.65)
ax_a.annotate("", xy=(0.89, 0.69), xytext=(0.89, 0.51),
              arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.1))
ax_a.annotate("", xy=(0.96, 0.51), xytext=(0.89, 0.51),
              arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.1))
ax_a.text(0.875, 0.72, r"$n=0$", color=BLUE, ha="center")
ax_a.text(0.95, 0.44, r"$n=-1$", color=ORANGE, ha="center")
ax_a.text(0.90, 0.25, r"$r_n(k_x,f)$", ha="center", fontweight="bold")

# b, angular resolution imposed by the beam waist.
delta_theta = np.linspace(-10, 10, 1001)
k0 = 2 * np.pi * F0 / C0
beam_colors = {25.: ORANGE, 100.: BLUE, 200.: PURPLE}
for waist in [25., 100., 200.]:
    delta_kx = k0 * (np.sin(np.deg2rad(THETA_BIC + delta_theta))
                     - np.sin(np.deg2rad(THETA_BIC)))
    intensity = np.exp(-0.5 * (waist * 1e-3 * delta_kx) ** 2)
    fwhm = np.rad2deg(GAUSSIAN_FWHM_FACTOR /
                     (k0 * waist * 1e-3 * np.cos(np.deg2rad(THETA_BIC))))
    label = rf"{int(waist)} mm  ({fwhm:.1f}$^\circ$)"
    ax_b.plot(delta_theta, intensity, color=beam_colors[waist], label=label)
ax_b.axvline(0, color=LIGHT, lw=0.7)
ax_b.set(xlabel=r"Angular offset, $\Delta\theta$ (deg)",
         ylabel="Incident intensity", xlim=(-10, 10), ylim=(0, 1.04),
         yticks=[0, 0.5, 1])
ax_b.legend(title="Waist (angular FWHM)", loc="upper right",
            handlelength=1.3, title_fontsize=7.2)
ax_b.text(-9.6, 0.87, "current", color=ORANGE, fontsize=7,
          fontweight="bold")

# c-d, ideal pole spectrum and the current finite-beam convolution.
extent = [theta.min(), theta.max(), frequency.min() / 1e3,
          frequency.max() / 1e3]
ideal_log = np.log10(np.maximum(ideal_map, 1e-6))
gaussian_log = np.log10(np.maximum(gaussian_maps[25.], 1e-6))
image = ax_c.imshow(ideal_log, origin="lower", extent=extent, aspect="auto",
                    cmap="magma", vmin=-6, vmax=0,
                    interpolation="bilinear", rasterized=True)
ax_d.imshow(gaussian_log, origin="lower", extent=extent, aspect="auto",
            cmap="magma", vmin=-6, vmax=0,
            interpolation="bilinear", rasterized=True)
for ax, title in [(ax_c, "Bloch plane wave"),
                  (ax_d, r"Gaussian beam, $w=25$ mm")]:
    ax.axvline(THETA_BIC, color="white", lw=0.75, ls=(0, (3, 2)))
    ax.set(xlabel=r"Central angle, $\theta$ (deg)", xlim=(5.8, 8.5),
           ylim=(179.55, 181.10))
    ax.set_title(title, pad=3)
ax_c.set_ylabel("Frequency (kHz)")
ax_d.tick_params(labelleft=False)
cbar = fig.colorbar(image, ax=[ax_c, ax_d], fraction=0.022, pad=0.012,
                    ticks=[-6, -3, 0])
cbar.set_label(r"Pole spectral density, $\log_{10}\mathcal{L}$")

# e, the intrinsic simultaneous closure of both homogeneous channels.
ax_e.semilogy(radiation["theta_deg"],
              np.maximum(radiation["radiation_A0"], 1e-17),
              color=BLUE, label=r"$|A_0|$")
ax_e.semilogy(radiation["theta_deg"],
              np.maximum(radiation["radiation_Am1"], 1e-17),
              color=ORANGE, label=r"$|A_{-1}|$")
ax_e.axvline(THETA_BIC, color=LIGHT, lw=0.7)
ax_e.set(xlabel=r"$\theta$ (deg)", ylabel="Eigenmode radiation",
         xlim=(7.099, 8.52), ylim=(1e-17, 1e-2),
         xticks=[7.1, 7.5, 8.0, 8.5])
ax_e.legend(loc="lower right", handlelength=1.3)

# f, intrinsic linewidth versus the angular-dispersion resolution floor.
ax_f.semilogy(theta, np.maximum(eigen_width, 1e-8), color=DARK,
              label="Bloch pole")
for waist in [25., 100., 200.]:
    ax_f.semilogy(theta, angular_bandwidths[waist],
                  color=beam_colors[waist], label=rf"$w={int(waist)}$ mm")
ax_f.axvline(THETA_BIC, color=LIGHT, lw=0.7)
ax_f.set(xlabel=r"Central angle, $\theta$ (deg)",
         ylabel="Spectral width (Hz)", xlim=(5.8, 8.5), ylim=(1e-8, 3e3))
ax_f.legend(loc="lower left", handlelength=1.2)

# g, experimentally actionable aperture/waist scaling.
minimum_width = np.array([np.nanmin(angular_bandwidths[w]) for w in waist_grid])
ax_g.loglog(waist_grid, minimum_width, "o-", color=BLUE,
            ms=3.5, mfc="white", mec=BLUE, mew=0.8)
ax_g.plot(25, minimum_width[0], "o", color=ORANGE, ms=5, zorder=3)
ax_g.annotate("current\n80-mm aperture", xy=(25, minimum_width[0]),
              xytext=(39, 780), color=ORANGE, fontsize=7,
              arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.7))
ax_g.set(xlabel="Beam waist (mm)", ylabel="Minimum width (Hz)",
         xlim=(20, 360), ylim=(80, 2e3))
ax_g.set_xticks([25, 50, 100, 200, 300])
ax_g.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_g.text(0.52, 0.12, r"matched aperture $\simeq3.2w$",
          transform=ax_g.transAxes, ha="center", va="bottom",
          color=GREY, fontsize=7)

# Shared styling and panel labels.
for label, ax in zip("abcdefg", [ax_a, ax_b, ax_c, ax_d, ax_e, ax_f, ax_g]):
    x_panel = -0.23 if ax in [ax_e, ax_f, ax_g] else -0.13
    y_panel = 1.15 if ax in [ax_e, ax_f, ax_g] else 1.08
    ax.text(x_panel, y_panel, label, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top")
for ax in [ax_b, ax_c, ax_d, ax_e, ax_f, ax_g]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

base = OUT / "Fig5_ideal_vs_finite_scattering"
fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
plt.close(fig)
