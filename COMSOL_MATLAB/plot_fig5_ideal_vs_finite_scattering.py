"""Figure 5: ideal Bloch scattering versus a finite Gaussian-beam sample."""
from pathlib import Path
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

ROOT = Path(__file__).resolve().parents[1]
FIN = ROOT / "COMSOL_MATLAB" / "rounded_180k_results" / "finite_20period_experiment"
INF = ROOT / "Ni2019_MATLAB" / "results" / "fig5_infinite_period_fields_180k"
POLE = ROOT / "Ni2019_MATLAB" / "results" / "fig4_wideangle_180k"
DRIVEN = ROOT / "Ni2019_MATLAB" / "results" / "fig4_experimental_observables_180k"
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


def continuous_spectrum(path, frequency_hz):
    """Dense finite-window Fourier integral; resolution remains ~2*pi/L."""
    data = read_csv(path)
    x = data["x_m"]
    field = data["reflected_real_pa"] + 1j * data["reflected_imag_pa"]
    k0 = 2 * np.pi * frequency_hz / 1500
    q = np.linspace(-1.25, 1.25, 1201)
    spectrum = np.exp(-1j * np.outer(q * k0, x)) @ (
        field * np.hanning(field.size))
    power = np.abs(spectrum) ** 2
    return q, 10 * np.log10(np.maximum(power / np.max(power), 1e-6))


def measured_fwhm(path):
    """FWHM of the background-subtracted ideal plane-wave Am1 response."""
    data = read_csv(path)
    angle_out, width_out = [], []
    for angle in np.unique(data["theta_deg"]):
        select = np.isclose(data["theta_deg"], angle, rtol=0, atol=2e-10)
        frequency = data["frequency_hz"][select]
        amplitude = (data["Am1_real"][select]
                     + 1j * data["Am1_imag"][select])
        valid = (np.isfinite(frequency) & np.isfinite(amplitude.real)
                 & np.isfinite(amplitude.imag))
        if np.count_nonzero(valid) < 15:
            continue
        frequency, amplitude = frequency[valid], amplitude[valid]
        order = np.argsort(frequency)
        frequency, amplitude = frequency[order], amplitude[order]
        intensity = np.abs(amplitude - 0.5 * (amplitude[0] + amplitude[-1])) ** 2
        peak = int(np.argmax(intensity))
        half = intensity[peak] / 2
        left = np.flatnonzero(intensity[:peak] <= half)
        right = np.flatnonzero(intensity[peak + 1:] <= half)
        if left.size == 0 or right.size == 0:
            continue
        il = left[-1]
        ir = peak + 1 + right[0]
        fl = np.interp(half, intensity[il:il + 2], frequency[il:il + 2])
        fr = np.interp(half, intensity[ir - 1:ir + 1][::-1],
                       frequency[ir - 1:ir + 1][::-1])
        if fr > fl:
            angle_out.append(angle)
            width_out.append(fr - fl)
    return np.array(angle_out), np.array(width_out)


def ideal_and_gaussian_maps(pole):
    """Ideal Lorentzian spectral density and its Gaussian angular average."""
    theta_positive = pole["theta_deg"]
    frequency_positive = pole["frequency_hz"]
    linewidth_positive = pole["linewidth_hz"]
    theta_internal = np.linspace(-20, 20, 2001)
    frequency_internal = np.interp(
        np.abs(theta_internal), theta_positive, frequency_positive)
    linewidth_internal = np.interp(
        np.abs(theta_internal), theta_positive, linewidth_positive)
    frequency_grid = np.linspace(179600, 180900, 1301)
    half = 0.5 * linewidth_internal
    ideal_internal = half[None, :] ** 2 / (
        (frequency_grid[:, None] - frequency_internal[None, :]) ** 2
        + half[None, :] ** 2 + 1e-300)

    theta_command = np.linspace(5.8, 8.5, 136)
    ideal_map = np.empty((frequency_grid.size, theta_command.size))
    for row in range(frequency_grid.size):
        ideal_map[row] = np.interp(
            theta_command, theta_internal, ideal_internal[row])

    # Source amplitude exp[-(x/w)^2] has angular intensity
    # exp[-(w*Delta kx)^2/2]. The 80-mm truncation is negligible at w=25 mm.
    waist = 25e-3
    k_reference = 2 * np.pi * 180e3 / 1500
    finite_map = np.empty_like(ideal_map)
    for column, theta_center in enumerate(theta_command):
        delta_kx = k_reference * (
            np.sin(np.deg2rad(theta_internal))
            - np.sin(np.deg2rad(theta_center)))
        weight = np.exp(-0.5 * (waist * delta_kx) ** 2)
        weight /= np.trapezoid(weight, theta_internal)
        finite_map[:, column] = np.trapezoid(
            ideal_internal * weight[None, :], theta_internal, axis=1)

    apparent_width = np.full(theta_command.size, np.nan)
    for column in range(theta_command.size):
        response = finite_map[:, column]
        above = np.flatnonzero(response >= 0.5 * np.max(response))
        if above.size > 1:
            apparent_width[column] = (
                frequency_grid[above[-1]] - frequency_grid[above[0]])
    return theta_command, frequency_grid, ideal_map, finite_map, apparent_width


finite_stems = ["gaussian_state_1_gamma", "gaussian_state_2_rayleigh_bic",
                "gaussian_state_3_above_bic"]
infinite_stems = ["gamma_limit", "rayleigh_bic_limit", "above_bic"]
frequencies = np.array([178780.2005, 180046.681327, 180498.322])
titles = [r"$0^\circ$", "Rayleigh BIC", r"$10.05^\circ$"]
colors = [GREY, BLUE, ORANGE]

xf = np.loadtxt(FIN / "finite_bic_field_x_m.csv", delimiter=",") * 1e3
yf = np.loadtxt(FIN / "finite_bic_field_y_m.csv", delimiter=",") * 1e3
xi = np.loadtxt(INF / "x_over_a.csv", delimiter=",")
yi = np.loadtxt(INF / "y_over_a.csv", delimiter=",")

finite_incident, finite_reflected, infinite_reflected = [], [], []
for finite_stem, infinite_stem in zip(finite_stems, infinite_stems):
    finite_incident.append(
        np.loadtxt(FIN / f"{finite_stem}_incident_real_pa.csv", delimiter=",")
        + 1j * np.loadtxt(FIN / f"{finite_stem}_incident_imag_pa.csv", delimiter=","))
    finite_reflected.append(
        np.loadtxt(FIN / f"{finite_stem}_reflected_real_pa.csv", delimiter=",")
        + 1j * np.loadtxt(FIN / f"{finite_stem}_reflected_imag_pa.csv", delimiter=","))
    infinite_reflected.append(
        np.loadtxt(INF / f"{infinite_stem}_reflected_real.csv", delimiter=",")
        + 1j * np.loadtxt(INF / f"{infinite_stem}_reflected_imag.csv", delimiter=","))

incident_scale = np.percentile(np.concatenate([
    np.abs(field[np.isfinite(field)]) for field in finite_incident]), 99.5)
finite_reflected = [field / incident_scale for field in finite_reflected]

orders = read_csv(INF / "floquet_orders.csv")
pole = read_csv(POLE / "wideangle_physical_pole.csv")
theta_scattering, width_scattering = measured_fwhm(
    DRIVEN / "adaptive_scattering_response.csv")
(theta_map, frequency_map, ideal_map, finite_map,
 finite_width) = ideal_and_gaussian_maps(pole)

# Figure contract: the ideal BIC is resolved by a Bloch plane wave and discrete
# Floquet projection, whereas a finite Gaussian beam angularly averages nearby
# leaky states and replaces the linewidth collapse by an aperture-limited floor.
fig = plt.figure(figsize=(7.2047244, 8.05))
outer = fig.add_gridspec(4, 1, height_ratios=[0.72, 0.72, 0.90, 0.92], hspace=0.52)
ideal_grid = outer[0].subgridspec(1, 3, wspace=0.08)
finite_grid = outer[1].subgridspec(1, 3, wspace=0.08)
ideal_axes = np.array([fig.add_subplot(ideal_grid[0, col]) for col in range(3)])
finite_axes = np.array([fig.add_subplot(finite_grid[0, col]) for col in range(3)])
middle = outer[2].subgridspec(1, 2, wspace=0.34)
ax_c = fig.add_subplot(middle[0, 0])
ax_d = fig.add_subplot(middle[0, 1])
bottom = outer[3].subgridspec(1, 3, wspace=0.40)
ax_e = fig.add_subplot(bottom[0, 0])
ax_f = fig.add_subplot(bottom[0, 1])
ax_g = fig.add_subplot(bottom[0, 2])

field_limit = 2.0
field_image = None
use_yi = (yi >= 0.12) & (yi <= 5.0)
extent_i = [xi.min(), xi.max(), yi[use_yi].min(), yi[use_yi].max()]
for column, ax in enumerate(ideal_axes):
    field_image = ax.imshow(np.real(infinite_reflected[column][use_yi]),
                            origin="lower", extent=extent_i, aspect="auto",
                            cmap="RdBu_r", vmin=-field_limit, vmax=field_limit,
                            interpolation="bilinear", rasterized=True)
    ax.axhline(0.12, color=DARK, lw=0.6)
    ax.set(xlim=(-2.5, 2.5), ylim=(0.12, 5), xticks=[-2, 0, 2], yticks=[1, 3, 5],
           xlabel="x/a")
    ax.set_title(titles[column], pad=3)
    if column:
        ax.tick_params(labelleft=False)
    else:
        ax.set_ylabel("y/a")

use_yf = (yf >= 0) & (yf <= 50)
extent_f = [xf.min(), xf.max(), 0, 50]
for column, ax in enumerate(finite_axes):
    ax.imshow(np.real(finite_reflected[column][use_yf]), origin="lower",
              extent=extent_f, aspect="auto", cmap="RdBu_r",
              vmin=-field_limit, vmax=field_limit,
              interpolation="bilinear", rasterized=True)
    ax.axhline(20, color=ORANGE, lw=0.6, ls=(0, (3, 2)))
    ax.set(xlim=(0, 148.4), ylim=(0, 50), xticks=[0, 50, 100, 150],
           yticks=[0, 20, 40], xlabel="x (mm)")
    if column:
        ax.tick_params(labelleft=False)
    else:
        ax.set_ylabel("y (mm)")
fig.colorbar(field_image, ax=[*ideal_axes, *finite_axes], fraction=0.010,
             pad=0.010, label=r"Re$(p_{\rm refl})/p_{\rm inc}$")

# c, exact discrete Floquet powers of the infinite periodic fields.
floor = 1e-8
for state, color, label in zip([1, 2, 3], colors, titles):
    select = (orders["state"] == state) & np.isin(orders["order"], [-1, 0])
    q = -orders["kx_over_k0"][select]
    eta = np.maximum(orders["power_fraction"][select], floor)
    visible = ((orders["order"][select] == 0)
               | (orders["is_propagating"][select] > 0))
    q, eta = q[visible], eta[visible]
    ax_c.vlines(q, floor, eta, color=color, lw=1.0, alpha=0.85)
    ax_c.plot(q, eta, "o", ms=4.0, mfc="white", mec=color, mew=0.9,
              label=label)
ax_c.set_yscale("log")
ax_c.axvline(-1, color=LIGHT, lw=0.65)
ax_c.axvline(1, color=LIGHT, lw=0.65)
ax_c.set(xlabel=r"Floquet line, $k_x/k_0$", ylabel="Power fraction",
         xlim=(-1.15, 1.15), ylim=(1e-8, 2), xticks=[-1, 0, 1])
ax_c.legend(loc="lower center", ncol=1, handlelength=0.8)

# d, continuous finite-aperture spectra extracted from the COMSOL scan line.
for stem, frequency, color, label in zip(
        finite_stems, frequencies, colors, titles):
    q, spectrum = continuous_spectrum(FIN / f"{stem}_scanline.csv", frequency)
    ax_d.plot(q, spectrum, color=color, label=label)
ax_d.axvspan(-1, 1, color="#F2F4F6", zorder=-3)
ax_d.axvline(-1, color=LIGHT, lw=0.65)
ax_d.axvline(1, color=LIGHT, lw=0.65)
ax_d.set(xlabel=r"Angular spectrum, $k_x/k_0$",
         ylabel="Windowed power (dB)", xlim=(-1.25, 1.25), ylim=(-60, 2))
ax_d.legend(loc="upper right", handlelength=1.4)

# e-f, ideal spectral pinch versus the finite-beam angular average.
ideal_log = np.log10(np.maximum(ideal_map, 1e-6))
finite_log = np.log10(np.maximum(finite_map, 1e-6))
extent_map = [theta_map.min(), theta_map.max(),
              frequency_map.min() / 1e3, frequency_map.max() / 1e3]
map_image = ax_e.imshow(ideal_log, origin="lower", extent=extent_map,
                        aspect="auto", cmap="magma", vmin=-6, vmax=0,
                        interpolation="bilinear", rasterized=True)
ax_f.imshow(finite_log, origin="lower", extent=extent_map, aspect="auto",
            cmap="magma", vmin=-6, vmax=0, interpolation="bilinear",
            rasterized=True)
for ax, title in [(ax_e, "Bloch plane wave"), (ax_f, "Gaussian beam")]:
    ax.axvline(7.0996629, color="white", lw=0.7, ls=(0, (3, 2)))
    ax.set(xlabel=r"$\theta$ (deg)", xlim=(5.8, 8.5), ylim=(179.6, 180.9))
    ax.set_title(title, pad=3)
ax_e.set_ylabel("f (kHz)")
ax_f.tick_params(labelleft=False)
cax = inset_axes(ax_f, width="38%", height="4%", loc="upper right", borderpad=0.8)
cbar = fig.colorbar(map_image, cax=cax, orientation="horizontal", ticks=[-6, 0])
cbar.ax.tick_params(labelsize=7.2, length=1.8, pad=1)
cbar.ax.set_title(r"$\log_{10}\mathcal{L}$", fontsize=7.2, pad=2, loc="left")

# g, field-extracted ideal linewidth and finite-beam apparent linewidth floor.
theta_bic = 7.0996629
theta_pole = pole["theta_deg"]
linewidth_pole = pole["linewidth_hz"]
show_pole = (theta_pole >= 5.8) & (theta_pole <= 8.5)
ax_g.semilogy(theta_pole[show_pole],
              np.maximum(linewidth_pole[show_pole], 1e-8),
              color=DARK, label="Eigenpole")
show_scattering = ((theta_scattering >= theta_bic + 0.015)
                   & (theta_scattering <= 8.5))
indices = np.flatnonzero(show_scattering)[::9]
ax_g.semilogy(theta_scattering[indices], width_scattering[indices], "o",
              ms=3.5, mfc="white", mec=ORANGE, mew=0.9,
              label="Plane wave")
ax_g.semilogy(theta_map, finite_width, color=BLUE, ls=(0, (3, 2)),
              label="Gaussian beam")
ax_g.axvline(theta_bic, color=LIGHT, lw=0.7, zorder=-1)
ax_g.set(xlabel=r"$\theta$ (deg)", ylabel="Extracted FWHM (Hz)",
         xlim=(5.8, 8.5), ylim=(1e-8, 2e3))
ax_g.legend(loc="lower left", handlelength=1.4)

# Panel labels and section headers.
ideal_axes[0].text(-0.20, 1.20, "a", transform=ideal_axes[0].transAxes,
                   fontsize=9, fontweight="bold", va="top")
ideal_axes[0].text(0.0, 1.20, "Ideal periodic scattering field",
                   transform=ideal_axes[0].transAxes, fontsize=8,
                   fontweight="bold", va="top")
finite_axes[0].text(-0.20, 1.20, "b", transform=finite_axes[0].transAxes,
                    fontsize=9, fontweight="bold", va="top")
finite_axes[0].text(0.0, 1.20, "Finite-beam scattering field",
                    transform=finite_axes[0].transAxes, fontsize=8,
                    fontweight="bold", va="top")
for label, ax in zip(["c", "d", "e", "f", "g"],
                     [ax_c, ax_d, ax_e, ax_f, ax_g]):
    ax.text(-0.18, 1.08, label, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top")
for ax in [*ideal_axes, *finite_axes, ax_c, ax_d, ax_e, ax_f, ax_g]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

base = OUT / "Fig5_ideal_vs_finite_scattering"
fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
plt.close(fig)
