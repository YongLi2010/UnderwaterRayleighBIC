"""Finite Gaussian-beam versus infinite Bloch-plane-wave scattering."""
from pathlib import Path
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIN = ROOT / "COMSOL_MATLAB" / "rounded_180k_results" / "finite_20period_experiment"
INF = ROOT / "Ni2019_MATLAB" / "results" / "fig4_infinite_period_fields_180k"
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
    data = read_csv(path)
    angle_out, width_out = [], []
    for angle in np.unique(data["theta_deg"]):
        sel = np.isclose(data["theta_deg"], angle, rtol=0, atol=2e-10)
        f = data["frequency_hz"][sel]
        a = data["Am1_real"][sel] + 1j * data["Am1_imag"][sel]
        valid = np.isfinite(f) & np.isfinite(a.real) & np.isfinite(a.imag)
        if np.count_nonzero(valid) < 15:
            continue
        f, a = f[valid], a[valid]
        order = np.argsort(f)
        f, a = f[order], a[order]
        response = np.abs(a - 0.5 * (a[0] + a[-1])) ** 2
        peak = int(np.argmax(response))
        half = response[peak] / 2
        left = np.flatnonzero(response[:peak] <= half)
        right = np.flatnonzero(response[peak + 1:] <= half)
        if left.size == 0 or right.size == 0:
            continue
        il = left[-1]
        ir = peak + 1 + right[0]
        fl = np.interp(half, response[il:il + 2], f[il:il + 2])
        fr = np.interp(half, response[ir - 1:ir + 1][::-1],
                       f[ir - 1:ir + 1][::-1])
        if fr > fl:
            angle_out.append(angle)
            width_out.append(fr - fl)
    return np.array(angle_out), np.array(width_out)


finite_stems = ["gaussian_state_1_gamma", "gaussian_state_2_rayleigh_bic",
                "gaussian_state_3_above_bic"]
infinite_stems = ["gamma_limit", "rayleigh_bic_limit", "above_bic"]
angles = np.array([0.0, 7.05, 10.05404])
frequencies = np.array([178780.2005, 180046.681327, 180498.322])
titles = [r"$0^\circ$", r"Rayleigh BIC", r"$10.05^\circ$"]
colors = [GREY, BLUE, ORANGE]

xf = np.loadtxt(FIN / "finite_bic_field_x_m.csv", delimiter=",") * 1e3
yf = np.loadtxt(FIN / "finite_bic_field_y_m.csv", delimiter=",") * 1e3
xi = np.loadtxt(INF / "x_over_a.csv", delimiter=",")
yi = np.loadtxt(INF / "y_over_a.csv", delimiter=",")

finite_incident, finite_reflected = [], []
infinite_incident, infinite_reflected = [], []
for fs, ins in zip(finite_stems, infinite_stems):
    finite_incident.append(
        np.loadtxt(FIN / f"{fs}_incident_real_pa.csv", delimiter=",")
        + 1j * np.loadtxt(FIN / f"{fs}_incident_imag_pa.csv", delimiter=","))
    finite_reflected.append(
        np.loadtxt(FIN / f"{fs}_reflected_real_pa.csv", delimiter=",")
        + 1j * np.loadtxt(FIN / f"{fs}_reflected_imag_pa.csv", delimiter=","))
    infinite_incident.append(
        np.loadtxt(INF / f"{ins}_incident_real.csv", delimiter=",")
        + 1j * np.loadtxt(INF / f"{ins}_incident_imag.csv", delimiter=","))
    infinite_reflected.append(
        np.loadtxt(INF / f"{ins}_reflected_real.csv", delimiter=",")
        + 1j * np.loadtxt(INF / f"{ins}_reflected_imag.csv", delimiter=","))

incident_scale = np.percentile(np.concatenate([
    np.abs(v[np.isfinite(v)]) for v in finite_incident]), 99.5)
finite_incident = [v / incident_scale for v in finite_incident]
finite_reflected = [v / incident_scale for v in finite_reflected]

orders = read_csv(INF / "floquet_orders.csv")
summary = read_csv(INF / "state_summary.csv")
pole = read_csv(POLE / "wideangle_physical_pole.csv")
theta_fwhm, width_fwhm = measured_fwhm(DRIVEN / "adaptive_scattering_response.csv")

# Figure contract: finite-aperture spectra are smooth but broadened, whereas
# an infinite Bloch plane wave produces discrete Floquet orders; neither fact
# alone proves a BIC, so the scattering-extracted linewidth is shown separately.
fig = plt.figure(figsize=(7.2047244, 7.45))
outer = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.72], hspace=0.35)
finite_grid = outer[0].subgridspec(2, 3, hspace=0.08, wspace=0.08)
infinite_grid = outer[1].subgridspec(2, 3, hspace=0.08, wspace=0.08)
fa = np.array([[fig.add_subplot(finite_grid[r, c]) for c in range(3)]
               for r in range(2)])
ia = np.array([[fig.add_subplot(infinite_grid[r, c]) for c in range(3)]
               for r in range(2)])
bottom = outer[2].subgridspec(1, 3, width_ratios=[1.25, 0.9, 0.95], wspace=0.40)
ax_c = fig.add_subplot(bottom[0, 0])
ax_d = fig.add_subplot(bottom[0, 1])
ax_e = fig.add_subplot(bottom[0, 2])

field_limit = 2.0
field_image = None
use_yf = (yf >= 0) & (yf <= 50)
extent_f = [xf.min(), xf.max(), 0, 50]
for col in range(3):
    for row, fields in enumerate([finite_incident, finite_reflected]):
        ax = fa[row, col]
        field_image = ax.imshow(np.real(fields[col][use_yf, :]), origin="lower",
                                extent=extent_f, aspect="auto", cmap="RdBu_r",
                                vmin=-field_limit, vmax=field_limit,
                                interpolation="bilinear", rasterized=True)
        ax.axhline(20, color=ORANGE, lw=0.6, ls=(0, (3, 2)))
        ax.set(xlim=(0, 148.4), ylim=(0, 50), xticks=[0, 50, 100, 150],
               yticks=[0, 20, 40])
        if row == 0:
            ax.set_title(titles[col], pad=3)
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("x (mm)")
        if col:
            ax.tick_params(labelleft=False)
        else:
            ax.set_ylabel("y (mm)")

use_yi = (yi >= 0.12) & (yi <= 5.0)
extent_i = [xi.min(), xi.max(), yi[use_yi].min(), yi[use_yi].max()]
for col in range(3):
    for row, fields in enumerate([infinite_incident, infinite_reflected]):
        ax = ia[row, col]
        ax.imshow(np.real(fields[col][use_yi, :]), origin="lower",
                  extent=extent_i, aspect="auto", cmap="RdBu_r",
                  vmin=-field_limit, vmax=field_limit,
                  interpolation="bilinear", rasterized=True)
        ax.axhline(0.12, color=DARK, lw=0.6)
        ax.set(xlim=(-2.5, 2.5), ylim=(0.12, 5), xticks=[-2, 0, 2],
               yticks=[1, 3, 5])
        if row == 0:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("x/a")
        if col:
            ax.tick_params(labelleft=False)
        else:
            ax.set_ylabel("y/a")

fa[0, 0].text(-0.30, 0.5, "Incident", rotation=90, rotation_mode="anchor",
              transform=fa[0, 0].transAxes, ha="center", va="center")
fa[1, 0].text(-0.30, 0.5, "Reflected", rotation=90, rotation_mode="anchor",
              transform=fa[1, 0].transAxes, ha="center", va="center")
ia[0, 0].text(-0.30, 0.5, "Incident", rotation=90, rotation_mode="anchor",
              transform=ia[0, 0].transAxes, ha="center", va="center")
ia[1, 0].text(-0.30, 0.5, "Reflected", rotation=90, rotation_mode="anchor",
              transform=ia[1, 0].transAxes, ha="center", va="center")
fig.colorbar(field_image, ax=[*fa.ravel(), *ia.ravel()], fraction=0.010, pad=0.010,
             label=r"Re$(p)/p_{\rm inc}$")

# c, continuous finite-window transform: smooth but resolution-limited.
for stem, freq, color, label in zip(finite_stems, frequencies, colors, titles):
    q, spectrum = continuous_spectrum(FIN / f"{stem}_scanline.csv", freq)
    ax_c.plot(q, spectrum, color=color, label=label)
ax_c.axvspan(-1, 1, color="#F2F4F6", zorder=-3)
ax_c.axvline(-1, color=LIGHT, lw=0.65)
ax_c.axvline(1, color=LIGHT, lw=0.65)
ax_c.set(xlabel=r"$k_x/k_0$", ylabel="Windowed power (dB)",
         xlim=(-1.25, 1.25), ylim=(-60, 2))
ax_c.legend(loc="upper right", ncol=1, handlelength=1.5)

# d, exact discrete far-field powers of the infinite periodic calculation.
floor = 1e-8
for state, color, label in zip([1, 2, 3], colors, titles):
    select = (orders["state"] == state) & np.isin(orders["order"], [-1, 0])
    q = -orders["kx_over_k0"][select]
    eta = np.maximum(orders["power_fraction"][select], floor)
    visible = (orders["order"][select] == 0) | (orders["is_propagating"][select] > 0)
    q, eta = q[visible], eta[visible]
    ax_d.vlines(q, floor, eta, color=color, lw=1.0, alpha=0.85)
    ax_d.plot(q, eta, "o", ms=4.0, mfc="white", mec=color, mew=0.9,
              label=label)
ax_d.set_yscale("log")
ax_d.axvline(-1, color=LIGHT, lw=0.65)
ax_d.axvline(1, color=LIGHT, lw=0.65)
ax_d.set(xlabel=r"Floquet line, $k_x/k_0$", ylabel="Power fraction",
         xlim=(-1.15, 1.15), ylim=(1e-8, 2), xticks=[-1, 0, 1])
ax_d.legend(loc="lower center", ncol=1, handlelength=0.8)

# e, the BIC discriminator: reflected-field FWHM matches the outgoing pole.
theta_bic = 7.0996629
theta = pole["theta_deg"]
linewidth = pole["linewidth_hz"]
open_side = (theta >= theta_bic) & (theta <= 8.55)
ax_e.semilogy(theta[open_side], np.maximum(linewidth[open_side], 1e-6),
              color=DARK, label="Eigenpole")
keep = (theta_fwhm >= theta_bic + 0.015) & (theta_fwhm <= 8.52)
indices = np.flatnonzero(keep)[::9]
ax_e.semilogy(theta_fwhm[indices], width_fwhm[indices], "o", ms=3.7,
              mfc="white", mec=ORANGE, mew=0.9, label="Scattering")
ax_e.axvline(theta_bic, color=LIGHT, lw=0.7, zorder=-1)
ax_e.set(xlabel=r"$\theta$ (deg)", ylabel="FWHM (Hz)",
         xlim=(7.08, 8.55), ylim=(1e-4, 8))
ax_e.legend(loc="lower right", handlelength=1.4)

fa[0, 0].text(-0.20, 1.28, "a", transform=fa[0, 0].transAxes,
              fontsize=9, fontweight="bold", va="top")
fa[0, 0].text(0.0, 1.28, "Finite Gaussian beam", transform=fa[0, 0].transAxes,
              fontsize=8, fontweight="bold", va="top")
ia[0, 0].text(-0.20, 1.28, "b", transform=ia[0, 0].transAxes,
              fontsize=9, fontweight="bold", va="top")
ia[0, 0].text(0.0, 1.28, "Infinite Bloch plane wave",
              transform=ia[0, 0].transAxes, fontsize=8, fontweight="bold", va="top")
for label, ax in zip(["c", "d", "e"], [ax_c, ax_d, ax_e]):
    ax.text(-0.20, 1.08, label, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top")
for ax in [*fa.ravel(), *ia.ravel(), ax_c, ax_d, ax_e]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

base = OUT / "Fig4_finite_vs_infinite_scattering"
fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
plt.close(fig)
