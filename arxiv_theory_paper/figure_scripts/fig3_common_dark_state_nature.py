#!/usr/bin/env python3
"""Nature-style Fig. 3: common-dark-state origin of the Rayleigh BIC."""
from pathlib import Path
import csv

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "Ni2019_MATLAB/results/fig3_radiation_contributions_180k/radiation_contributions.csv"
OUT = ROOT / "arxiv_theory_paper/figures"

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.5, "axes.titlesize": 8.2, "axes.titleweight": "semibold",
    "axes.linewidth": 0.7, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.size": 2.6, "ytick.major.size": 2.6,
    "xtick.major.width": 0.65, "ytick.major.width": 0.65,
    "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
})

INK = "#202A32"; BLUE = "#21618C"; TEAL = "#168A8A"
GOLD = "#D88418"; RED = "#C9473B"; PURPLE = "#7263A5"
GREY = "#78858B"; PALE = "#F2F5F5"


def table(path):
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def z(d, name):
    return d[name + "_real"] + 1j*d[name + "_imag"]


def label(ax, s, x=-0.11, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="bottom", ha="left", color=INK, clip_on=False)


def bare(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)


def cmt_schematic(ax):
    bare(ax); label(ax, "a")
    ax.set_title("Two modes, two radiation channels", pad=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    # Reduced two-mode system, deliberately abstracted from the geometry.
    for x, col, txt in [(0.28, TEAL, "L"), (0.72, GOLD, "S")]:
        ax.add_patch(Circle((x, .45), .095, fc=col, ec="white", lw=.8))
        ax.text(x, .45, txt, color="white", ha="center", va="center",
                fontsize=10, fontweight="bold")
    ax.text(.28, .31, "large-cavity mode", ha="center", va="top",
            color=TEAL, fontsize=7.0)
    ax.text(.72, .31, "small-cavity mode", ha="center", va="top",
            color=GOLD, fontsize=7.0)
    ax.add_patch(FancyArrowPatch((.38,.45),(.62,.45),arrowstyle="<->",
                                 mutation_scale=9,lw=1.15,color=PURPLE))
    ax.text(.50,.52,r"$H(k_x)$",ha="center",color=PURPLE)
    for y, n, col in [(.75,"0",BLUE),(.13,"-1",RED)]:
        ax.add_patch(FancyBboxPatch((.40,y-.055),.20,.11,
                     boxstyle="round,pad=.012,rounding_size=.025",
                     fc=PALE,ec=col,lw=.9))
        ax.text(.50,y,rf"$J_{{{n}}}$",ha="center",va="center",color=col)
        for x0 in (.28,.72):
            ax.add_patch(FancyArrowPatch((x0,.55 if y>.5 else .35),(.45 if x0<.5 else .55,y-.06 if y>.5 else y+.06),
                         arrowstyle="-|>",mutation_scale=7,lw=.8,color=col,
                         connectionstyle="arc3,rad=0"))
    ax.text(.50,.91,r"$\mathbf{J}=D_s(k_x)\mathbf{u}$",ha="center",va="top",
            color=INK,fontsize=8.5)
    ax.text(.50,-.01,"one Bloch eigenvector must darken both rows",ha="center",
            color=GREY,fontsize=7.2)


def kstyle(ax, kb):
    ax.axvline(kb,color=RED,lw=.9,ls=(0,(2.5,2.0)),zorder=0)
    ax.set_xlim(.088,.132); ax.set_xticks([.09,.10,.11,.12,.13])
    ax.tick_params(labelsize=7.3,pad=1.8)
    ax.set_xlabel(r"Bloch wavevector $k_xa/2\pi$")


def hybrid(ax, d, kb, ib):
    label(ax,"b"); ax.set_title("Bloch phase tunes the hybrid coordinate",pad=5)
    al, ass = z(d,"a_large"), z(d,"a_small")
    norm = np.abs(al)**2 + np.abs(ass)**2
    wl, ws = np.abs(al)**2/norm, np.abs(ass)**2/norm
    ax.plot(d["kappa"],wl,color=TEAL,lw=1.5,label="Large mode")
    ax.plot(d["kappa"],ws,color=GOLD,lw=1.5,label="Small mode")
    ax.scatter([kb,kb],[wl[ib],ws[ib]],s=25,c=[TEAL,GOLD],ec="white",lw=.5,zorder=4)
    ax.set_ylim(.36,.64); ax.set_yticks([.4,.5,.6]); ax.set_ylabel("Modal weight")
    kstyle(ax,kb); ax.legend(loc="upper left",ncol=2,handlelength=1.4,columnspacing=1.0)


def collapse(ax, d, kb, ib):
    label(ax,"c"); ax.set_title("Eigenstate enters the common dark subspace",pad=5)
    floor=1e-17
    s=np.maximum(d["radiation_sigma_ratio"],floor)
    m=np.maximum(d["dark_state_misalignment"],floor)
    ax.semilogy(d["kappa"],s,color=BLUE,lw=1.55,label=r"$\sigma_{\min}(D_s)/\sigma_{\max}(D_s)$")
    ax.semilogy(d["kappa"],m,color=PURPLE,lw=1.45,label=r"$1-|\langle u_{\rm dark}|u\rangle|^2$")
    ax.scatter([kb,kb],[s[ib],m[ib]],s=25,c=[BLUE,PURPLE],ec="white",lw=.45,zorder=4)
    ax.set_ylim(1e-17,2); ax.set_yticks([1e-16,1e-12,1e-8,1e-4,1])
    ax.set_ylabel("Normalized residual"); kstyle(ax,kb)
    ax.legend(loc="upper center",fontsize=7.2,handlelength=1.5)


def phasor(ax, d, ib):
    bare(ax); label(ax,"d")
    ax.set_title("Simultaneous channel cancellation",pad=5)
    ax.set_xlim(-1.25,1.25); ax.set_ylim(-1.02,1.12); ax.set_aspect("equal",adjustable="box")
    pairs=[(z(d,"b0_large")[ib],z(d,"b0_small")[ib],.50,BLUE,r"$n=0$"),
           (z(d,"bm1_large")[ib],z(d,"bm1_small")[ib],-.50,RED,r"$n=-1$")]
    for bl,bs,y0,col,name in pairs:
        rot=np.exp(-1j*np.angle(bl)); scale=max(abs(bl),abs(bs))
        vl,vs=bl*rot/scale,bs*rot/scale
        ax.plot([-1.13,1.13],[y0,y0],color="#D6DCDE",lw=.65)
        ax.add_patch(FancyArrowPatch((0,y0),(vl.real,y0+vl.imag),arrowstyle="-|>",
                                     mutation_scale=9,lw=1.65,color=TEAL))
        ax.add_patch(FancyArrowPatch((0,y0),(vs.real,y0+vs.imag),arrowstyle="-|>",
                                     mutation_scale=9,lw=1.65,color=GOLD))
        end=vl+vs
        ax.scatter([end.real],[y0+end.imag],s=34,marker="x",color=col,lw=1.8,zorder=5)
        ax.text(-1.20,y0+.18,name,color=col,ha="left",fontweight="semibold")
    ax.text(.46,.82,"L",color=TEAL,fontweight="semibold")
    ax.text(-.56,.82,"S",color=GOLD,fontweight="semibold")
    ax.text(0,-.93,r"$J_n=J_{nL}+J_{nS}=0$",ha="center",color=INK,fontsize=8.2)


def main():
    d=table(SOURCE); k=d["kappa"]
    score=d["b0_sum_abs"]+d["bm1_sum_abs"]
    ib=int(np.argmin(score)); kb=k[ib]
    assert d["radiation_sigma_ratio"][ib] < 1e-12
    assert d["dark_state_misalignment"][ib] < 1e-12
    fig=plt.figure(figsize=(7.15,4.72),facecolor="white")
    gs=fig.add_gridspec(2,2,left=.07,right=.985,bottom=.11,top=.94,
                        wspace=.28,hspace=.38)
    cmt_schematic(fig.add_subplot(gs[0,0]))
    hybrid(fig.add_subplot(gs[0,1]),d,kb,ib)
    collapse(fig.add_subplot(gs[1,0]),d,kb,ib)
    phasor(fig.add_subplot(gs[1,1]),d,ib)
    OUT.mkdir(parents=True,exist_ok=True)
    stem=OUT/"fig3_common_dark_state_nature"
    fig.savefig(stem.with_suffix(".pdf"),bbox_inches="tight",facecolor="white")
    fig.savefig(stem.with_suffix(".svg"),bbox_inches="tight",facecolor="white")
    fig.savefig(stem.with_suffix(".png"),dpi=600,bbox_inches="tight",facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"),dpi=600,bbox_inches="tight",facecolor="white",
                pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig)
    print(stem)
    print(f"kappa_BIC={kb:.15g}; sigma ratio={d['radiation_sigma_ratio'][ib]:.3e}; misalignment={d['dark_state_misalignment'][ib]:.3e}")


if __name__ == "__main__": main()
