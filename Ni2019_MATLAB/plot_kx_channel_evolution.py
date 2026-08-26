#!/usr/bin/env python3
"""PRL-style radiation-channel evolution along the Rayleigh-BIC pole branch."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "kx_channel_evolution_180k"
CSV = OUT / "kx_channel_evolution.csv"

data = np.genfromtxt(CSV, delimiter=",", names=True)
kappa = data["kappa"]
dk = data["delta_kappa"]
omega = data["Omega_real"]
a0 = data["A0_over_C"]
am = data["Am1_over_C"]
k_bic = float(kappa[dk == 0][0])
omega_bic = float(omega[dk == 0][0])

fit = (np.abs(dk) >= 1e-6) & (np.abs(dk) <= 1e-3)
p0 = np.polyfit(np.log10(np.abs(dk[fit])), np.log10(a0[fit]), 1)
pm = np.polyfit(np.log10(np.abs(dk[fit])), np.log10(am[fit]), 1)

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 8.0,
        "axes.labelsize": 8.0,
        "xtick.labelsize": 7.4,
        "ytick.labelsize": 7.4,
        "axes.linewidth": 0.75,
        "xtick.major.width": 0.65,
        "ytick.major.width": 0.65,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "legend.fontsize": 7.2,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "savefig.transparent": False,
    }
)

navy = "#183B66"
orange = "#D97721"
red = "#C7352C"
gray = "#6F7478"
ev_fill = "#EAF1F7"
op_fill = "#FBF0E6"

fig, axes = plt.subplots(1, 3, figsize=(183 / 25.4, 70 / 25.4))
fig.subplots_adjust(left=0.063, right=0.992, bottom=0.205, top=0.905, wspace=0.36)


def shade_sides(ax):
    ax.axvspan(kappa.min(), k_bic, color=ev_fill, zorder=-5)
    ax.axvspan(k_bic, kappa.max(), color=op_fill, zorder=-5)
    ax.axvline(k_bic, color=red, lw=0.85, ls=(0, (2.2, 2.2)), zorder=1)


# a, pole branch crosses the n=-1 Rayleigh line.
ax = axes[0]
shade_sides(ax)
ax.plot(kappa, omega, color=navy, lw=1.65, label="Pole")
ax.plot(kappa, 1 - kappa, color=gray, lw=1.15, ls=(0, (3.2, 2.2)), label="Rayleigh line")
ax.scatter([k_bic], [omega_bic], s=22, color=red, edgecolor="white", linewidth=0.6, zorder=5)
ax.text(0.18, 0.93, "evanescent", color=navy, transform=ax.transAxes, ha="center", va="top", fontsize=7.2)
ax.text(0.82, 0.93, "open", color=orange, transform=ax.transAxes, ha="center", va="top", fontsize=7.2)
ax.set_xlim(kappa.min(), kappa.max())
ax.set_ylim(0.858, 0.920)
ax.set_xlabel(r"Bloch wave number, $\kappa=k_xa/2\pi$")
ax.set_ylabel(r"Normalized frequency, $\mathrm{Re}\,\Omega$")
ax.legend(loc="lower left", handlelength=2.2, borderaxespad=0.2)

# b, the two radiation amplitudes share one zero.
ax = axes[1]
shade_sides(ax)
left = dk < 0
right = dk > 0
for mask in (left, right):
    ax.semilogy(kappa[mask], a0[mask], color=navy, lw=1.65)
    ax.semilogy(kappa[mask], am[mask], color=orange, lw=1.65)
ax.annotate(
    r"$A_0=A_{-1}=0$",
    xy=(k_bic, 1.1e-7),
    xytext=(k_bic + 0.004, 3.0e-7),
    color=red,
    fontsize=7.4,
    ha="left",
    va="bottom",
    arrowprops=dict(arrowstyle="-|>", color=red, lw=0.75, shrinkA=1.5, shrinkB=1.5),
)
ax.text(0.96, 0.79, r"$|A_{-1}|/\|C\|_2$", color=orange, transform=ax.transAxes, ha="right", va="center")
ax.text(0.96, 0.55, r"$|A_0|/\|C\|_2$", color=navy, transform=ax.transAxes, ha="right", va="center")
ax.set_xlim(kappa.min(), kappa.max())
ax.set_ylim(8e-8, 2e-2)
ax.set_xlabel(r"Bloch wave number, $\kappa$")
ax.set_ylabel("Normalized radiation amplitude")

# c, common first-order opening of both channels.
ax = axes[2]
x = np.abs(dk[dk != 0])
for values, color, marker in ((a0, navy, "o"), (am, orange, "s")):
    yneg = values[left]
    ypos = values[right]
    xneg = np.abs(dk[left])
    xpos = np.abs(dk[right])
    ax.loglog(xneg, yneg, marker=marker, ms=2.7, mfc=color, mec=color, mew=0.4, ls="none", alpha=0.60)
    ax.loglog(xpos, ypos, marker=marker, ms=2.7, mfc="white", mec=color, mew=0.75, ls="none", alpha=0.95)

xfit = np.logspace(-6, -3, 80)
ax.loglog(xfit, 10 ** p0[1] * xfit ** p0[0], color=navy, lw=1.5)
ax.loglog(xfit, 10 ** pm[1] * xfit ** pm[0], color=orange, lw=1.5)
ax.text(0.05, 0.90, rf"$A_{{-1}}$: slope {pm[0]:.3f}", color=orange, transform=ax.transAxes, ha="left", va="top")
ax.text(0.05, 0.78, rf"$A_0$: slope {p0[0]:.3f}", color=navy, transform=ax.transAxes, ha="left", va="top")
ax.set_xlim(8e-7, 3.2e-2)
ax.set_ylim(3e-8, 2e-2)
ax.set_xlabel(r"Distance from BIC, $|\Delta\kappa|$")
ax.set_ylabel("Normalized radiation amplitude")

for label, ax in zip("abc", axes):
    ax.text(-0.18, 1.06, label, transform=ax.transAxes, fontsize=8.6, fontweight="bold", va="top")
    ax.tick_params(direction="out", pad=2.0)

stem = OUT / "kx_channel_evolution_prl"
fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"A0 exponent: {p0[0]:.8f}")
print(f"A-1 exponent: {pm[0]:.8f}")
print(stem)
