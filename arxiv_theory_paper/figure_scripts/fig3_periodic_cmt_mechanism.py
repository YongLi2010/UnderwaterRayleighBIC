#!/usr/bin/env python3
"""Fig. 3 | Periodic coupled-mode origin of the Rayleigh BIC.

The full modal-matching eigenbranch is projected onto two fixed dressed
aperture modes extracted at the strict BIC.  The wide-mode basis retains all
wide-groove transverse orders; likewise for the narrow mode.  The resulting
two-mode coordinates and the exact groove-resolved channel source terms are
then used to display simultaneous two-channel destructive interference.
"""
from pathlib import Path
import csv

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (ROOT / "Ni2019_MATLAB" / "results" /
          "fig3_radiation_contributions_180k" /
          "radiation_contributions.csv")
OUT = ROOT / "arxiv_theory_paper" / "figures"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.4,
    "axes.titlesize": 8.2,
    "axes.titleweight": "semibold",
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

NAVY = "#173A55"
TEAL = "#148C8C"
GOLD = "#D8891C"
PURPLE = "#7868A6"
RED = "#D94B3D"
GREY = "#6C757A"
PALE = "#EEF2F3"
LIGHT_TEAL = "#D9EEEE"
LIGHT_GOLD = "#F7E9D2"


def read_table(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    return {key: np.array([float(row[key]) for row in rows])
            for key in rows[0]}


def z(data, prefix):
    return data[f"{prefix}_real"] + 1j*data[f"{prefix}_imag"]


def panel_label(ax, text, x=-0.07, y=1.045):
    ax.text(x, y, text, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9.4, fontweight="bold", color="#111111")


def clear_axis(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_periodic_cmt(ax):
    clear_axis(ax); panel_label(ax, "a")
    ax.set_title("Periodic two-mode lattice", pad=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # Three unit cells of a diatomic resonator chain.  The central cell is the
    # visual anchor, while the faded neighbors encode periodic continuation.
    cells = [(0.15, 0.34), (0.42, 0.61), (0.69, 0.88)]
    for idx, (xl, xr) in enumerate(cells):
        alpha = 1.0 if idx == 1 else 0.45
        if idx == 1:
            ax.add_patch(FancyBboxPatch((xl-0.035, 0.28), xr-xl+0.07, 0.42,
                         boxstyle="round,pad=0.012,rounding_size=0.018",
                         facecolor="#F8FAFA", edgecolor="#CBD3D6", lw=0.8))
        xL, xS, y = xl+0.045, xr-0.045, 0.49
        ax.plot([xL+0.055, xS-0.055], [y, y], color=NAVY, lw=1.15,
                alpha=alpha)
        ax.add_patch(Circle((xL, y), 0.055, facecolor=TEAL,
                            edgecolor="white", lw=0.7, alpha=alpha))
        ax.add_patch(Circle((xS, y), 0.055, facecolor=GOLD,
                            edgecolor="white", lw=0.7, alpha=alpha))
        if idx == 1:
            ax.text(xL, y, "L", color="white", ha="center", va="center",
                    fontweight="bold")
            ax.text(xS, y, "S", color="white", ha="center", va="center",
                    fontweight="bold")
            ax.text((xL+xS)/2, y+0.055, r"$J_0$", ha="center", color=NAVY)
            ax.text((xl+xr)/2, 0.31, "cell $m$", ha="center", color=GREY)

    # Inter-cell couplings and periodic continuation.
    for x1, x2 in [(0.295, 0.465), (0.565, 0.735)]:
        ax.plot([x1, x2], [0.49, 0.49], color=PURPLE, lw=1.05, alpha=0.8)
        ax.text((x1+x2)/2, 0.545, r"$J_1$", ha="center", color=PURPLE)
    # Draw continuation dots as vector primitives so the PDF has no
    # font-dependent ellipsis glyph and remains visually symmetric.
    for xc in (0.020, 0.965):
        for dx in (-0.018, 0.0, 0.018):
            ax.add_patch(Circle((xc + dx, 0.49), 0.0065,
                                facecolor=GREY, edgecolor="none"))

    ax.annotate("", xy=(0.80, 0.83), xytext=(0.26, 0.83),
                arrowprops=dict(arrowstyle="-|>", color=NAVY, lw=1.05))
    ax.text(0.53, 0.89, r"Bloch phase $\exp(ik_xa)$", ha="center", color=NAVY)
    ax.plot([0.39, 0.64], [0.19, 0.19], color=RED, lw=0.9)
    ax.plot([0.39, 0.39], [0.17, 0.21], color=RED, lw=0.9)
    ax.plot([0.64, 0.64], [0.17, 0.21], color=RED, lw=0.9)
    ax.text(0.515, 0.12, "period $a$", ha="center", color=RED)
    ax.text(0.20, 0.73, "dressed cavity modes", ha="center", color=GREY)


def draw_bloch_map(ax):
    clear_axis(ax); panel_label(ax, "b")
    ax.set_title("Bloch-reduced coupled-mode theory", pad=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # kx controls two inseparable pieces of the reduced model: the hybrid
    # eigenvector and the pressure-channel dark-state requirements.
    ax.add_patch(Circle((0.11, 0.50), 0.065, facecolor=NAVY,
                        edgecolor="white", lw=0.7))
    ax.text(0.11, 0.50, r"$k_x$", color="white", ha="center", va="center",
            fontsize=8.2, fontweight="semibold")

    boxes = [
        (0.28, 0.64, 0.25, 0.19, LIGHT_TEAL, TEAL,
         r"$H_{\rm B}(k_x)$", "hybrid eigenstate"),
        (0.28, 0.17, 0.25, 0.19, LIGHT_GOLD, GOLD,
         r"$D_p(k_x)$", "channel dark states"),
    ]
    for x0, y0, w, h, fc, ec, math, caption in boxes:
        ax.add_patch(FancyBboxPatch((x0, y0), w, h,
                     boxstyle="round,pad=0.012,rounding_size=0.018",
                     facecolor=fc, edgecolor=ec, lw=0.9))
        ax.text(x0+w/2, y0+h*0.62, math, ha="center", va="center",
                color=ec, fontsize=8.2)
        ax.text(x0+w/2, y0+h*0.25, caption, ha="center", va="center",
                color=NAVY, fontsize=6.7)
        ax.add_patch(FancyArrowPatch((0.18, 0.50), (x0, y0+h/2),
                     connectionstyle="arc3,rad=0.0", arrowstyle="-|>",
                     mutation_scale=7.5, color="#8C969B", lw=0.9))

    outputs = [
        (0.76, 0.735, NAVY, r"$r_{\rm eig}$"),
        (0.76, 0.315, TEAL, r"$r_0$"),
        (0.76, 0.195, GOLD, r"$r_{-1}$"),
    ]
    for x0, y0, color, text in outputs:
        ax.add_patch(Circle((x0, y0), 0.054, facecolor="white",
                            edgecolor=color, lw=1.2))
        ax.text(x0, y0, text, ha="center", va="center", color=color,
                fontsize=7.5)
    ax.add_patch(FancyArrowPatch((0.53, 0.735), (0.70, 0.735),
                                 arrowstyle="-|>", mutation_scale=7.5,
                                 color=NAVY, lw=0.9))
    for y0 in (0.315, 0.195):
        ax.add_patch(FancyArrowPatch((0.53, 0.265), (0.70, y0),
                         connectionstyle="arc3,rad=0.0", arrowstyle="-|>",
                         mutation_scale=7.5, color="#8C969B", lw=0.85))

    # Visual common-null-state condition; the full equation belongs in text.
    ax.add_patch(FancyBboxPatch((0.84, 0.37), 0.13, 0.26,
                 boxstyle="round,pad=0.012,rounding_size=0.018",
                 facecolor="#FCEAE7", edgecolor=RED, lw=0.9))
    for y0, color in [(0.56, NAVY), (0.50, TEAL), (0.44, GOLD)]:
        ax.add_patch(FancyArrowPatch((0.80, y0), (0.875, 0.50),
                     arrowstyle="-|>", mutation_scale=6.8,
                     color=color, lw=0.85))
    ax.add_patch(Circle((0.90, 0.50), 0.020, facecolor=RED,
                        edgecolor="white", lw=0.4))
    ax.text(0.905, 0.31, "common dark state", ha="center", color=RED,
            fontsize=6.7, fontweight="semibold")


def style_k_axis(ax, k_bic, show_xlabel=True):
    ax.axvline(k_bic, color=RED, lw=0.85, ls=(0, (2.6, 2.2)), zorder=1)
    ax.tick_params(labelsize=6.6, pad=1.8)
    ax.set_xlim(0.088, 0.132)
    ax.set_xticks([0.09, 0.10, 0.11, 0.12, 0.13])
    if show_xlabel:
        ax.set_xlabel(r"Bloch wavevector $k_xa/2\pi$")


def draw_hybridization(fig, spec, k, weight_l, phase, k_bic):
    outer = fig.add_subplot(spec); clear_axis(outer); panel_label(outer, "b")
    outer.set_title("Bloch hybridization", pad=5)
    sub = spec.subgridspec(2, 1, hspace=0.12)
    aw = fig.add_subplot(sub[0]); ap = fig.add_subplot(sub[1], sharex=aw)
    aw.plot(k, weight_l, color=TEAL, lw=1.45)
    aw.plot(k, 1-weight_l, color=GOLD, lw=1.45)
    aw.set_ylabel("Modal\nweight", rotation=0, rotation_mode="anchor",
                  ha="right", va="center", labelpad=12, color=NAVY)
    aw.set_ylim(0.40, 0.59); aw.set_yticks([0.4, 0.5, 0.6])
    aw.text(k[-1]-0.001, weight_l[-1]-0.008, "L", color=TEAL,
            ha="right", va="top", fontweight="semibold")
    aw.text(k[-1]-0.001, 1-weight_l[-1]+0.007, "S", color=GOLD,
            ha="right", va="bottom", fontweight="semibold")
    ap.plot(k, phase, color=PURPLE, lw=1.45)
    ap.axhline(0, color="#C7CED1", lw=0.55)
    ap.set_ylabel(r"$\Delta\phi/\pi$", rotation=0, rotation_mode="anchor",
                  ha="right", va="center", labelpad=12, color=NAVY)
    lim=max(abs(phase))*1.18; ap.set_ylim(-lim,lim)
    for ax in (aw, ap): style_k_axis(ax, k_bic, ax is ap)
    aw.tick_params(labelbottom=False)
    ib=int(np.argmin(abs(k-k_bic)))
    aw.scatter([k_bic,k_bic],[weight_l[ib],1-weight_l[ib]],s=20,
               c=[TEAL,GOLD],edgecolors="white",linewidths=0.45,zorder=5)
    ap.scatter([k_bic],[phase[ib]],s=22,c=RED,edgecolors="white",
               linewidths=0.45,zorder=5)
    return outer


def matching_panel(ax, k, curves, k_bic, ylabel, title, letter,
                   annotation, phase=False):
    panel_label(ax, letter); ax.set_title(title, pad=5)
    styles = [(TEAL, "-", 1.55),
              (GOLD, (0, (2.0, 1.5)), 1.65)]
    for curve, (color, ls, lw) in zip(curves, styles):
        ax.plot(k, curve, color=color, ls=ls, lw=lw)
    style_k_axis(ax, k_bic, True)
    ax.axhline(1.0, color="#BFC7CA", lw=0.65, zorder=0)
    ax.set_ylabel(ylabel, color=NAVY)
    ib=int(np.argmin(abs(k-k_bic)))
    y_bic=np.mean([curve[ib] for curve in curves])
    ax.scatter([k_bic],[y_bic],s=31,c=RED,edgecolors="white",
               linewidths=0.55,zorder=6)
    ax.annotate(annotation, xy=(k_bic,y_bic),
                xytext=(k_bic+0.0032,
                        y_bic+(0.0045 if not phase else 0.0010)),
                color=RED, fontsize=6.8, ha="left", va="bottom",
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.65))
    ax.grid(False)


def main():
    data = read_table(SOURCE)
    k = data["kappa"]
    a_l, a_s = z(data, "a_large"), z(data, "a_small")
    b0_l, b0_s = z(data, "b0_large"), z(data, "b0_small")
    bm_l, bm_s = z(data, "bm1_large"), z(data, "bm1_small")
    closure = data["b0_sum_abs"] + data["bm1_sum_abs"]
    ib = int(np.argmin(closure)); k_bic = k[ib]

    weight_l = abs(a_l)**2/(abs(a_l)**2+abs(a_s)**2)
    hybrid_phase = np.unwrap(np.angle(a_s/a_l))/np.pi
    amplitude = [abs(b0_l)/abs(b0_s),abs(bm_l)/abs(bm_s)]
    phase=[]
    for ratio in (b0_l/b0_s,bm_l/bm_s):
        value=np.angle(ratio)/np.pi
        phase.append(np.where(value<0,value+2,value))

    # The fixed dressed basis must remain a faithful local two-mode subspace.
    assert np.min(data["basis_fidelity_large"]) > 0.995
    assert np.min(data["basis_fidelity_small"]) > 0.995
    assert closure[ib] < 1e-12

    fig = plt.figure(figsize=(7.0,5.15),facecolor="white")
    gs = fig.add_gridspec(2,2,left=0.065,right=0.985,bottom=0.095,
                          top=0.95,wspace=0.25,hspace=0.34,
                          height_ratios=[0.88,1.12])
    axa=fig.add_subplot(gs[0,0]); draw_periodic_cmt(axa)
    draw_hybridization(fig,gs[0,1],k,weight_l,hybrid_phase,k_bic)
    axc=fig.add_subplot(gs[1,0])
    matching_panel(axc,k,amplitude,k_bic,
                   r"Source ratio $|J_{nL}|/|J_{nS}|$",
                   "Amplitude balance","c","equal amplitudes",phase=False)
    axc.set_ylim(0.985,1.047); axc.set_yticks([0.99,1.00,1.02,1.04])
    axd=fig.add_subplot(gs[1,1])
    matching_panel(axd,k,phase,k_bic,r"Phase difference $\Delta\phi_n/\pi$",
                   "Phase opposition","d",r"$\pi$ phase shift",phase=True)
    axd.set_ylim(0.991,1.004); axd.set_yticks([0.992,0.996,1.000,1.004])

    handles=[Line2D([0],[0],color=TEAL,lw=1.55,label=r"$n=0$"),
             Line2D([0],[0],color=GOLD,lw=1.65,ls=(0,(2,1.5)),label=r"$n=-1$")]
    fig.legend(handles=handles,loc="center",bbox_to_anchor=(0.51,0.483),
               ncol=2,handlelength=2.2,columnspacing=1.7,fontsize=6.9)

    OUT.mkdir(parents=True,exist_ok=True)
    stem=OUT/"fig3_periodic_cmt_mechanism"
    fig.savefig(stem.with_suffix(".pdf"),bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"),bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"),dpi=600,bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"),dpi=600,bbox_inches="tight")
    plt.close(fig)
    print(stem)
    print(f"kappa_BIC={k_bic:.15g}, source closure={closure[ib]:.3e}")
    print("basis fidelity minima:",data["basis_fidelity_large"].min(),
          data["basis_fidelity_small"].min())


if __name__ == "__main__":
    main()
