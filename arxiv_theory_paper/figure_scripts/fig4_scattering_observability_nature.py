#!/usr/bin/env python3
"""Nature-style Fig. 4: experimentally accessible scattering signatures."""
from pathlib import Path
import csv

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle, Arc
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "Ni2019_MATLAB/results/fig4_experimental_observables_180k"
OUT = ROOT / "arxiv_theory_paper/figures"

mpl.rcParams.update({
    "font.family":"sans-serif", "font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "font.size":7.3, "axes.titlesize":8.1, "axes.titleweight":"semibold",
    "axes.linewidth":.7, "axes.spines.top":False, "axes.spines.right":False,
    "xtick.major.size":2.5, "ytick.major.size":2.5,
    "xtick.major.width":.65, "ytick.major.width":.65,
    "legend.frameon":False, "pdf.fonttype":42, "svg.fonttype":"none",
})

INK="#202A32"; BLUE="#245F87"; TEAL="#158B8A"; GOLD="#D88418"
RED="#C9473B"; PURPLE="#7162A3"; GREY="#7C898E"; WATER="#EDF5F5"
PHASE=LinearSegmentedColormap.from_list("phase",["#214E73","#94BED0","#F7F7F3","#EBA26A","#9F3A35"])
FIELD=LinearSegmentedColormap.from_list("field",["#214E73","#9BC5D3","#FAF8F1","#ECA06A","#A63D36"])
LINE=LinearSegmentedColormap.from_list("line",["#F7F7F3","#BBD8D6","#2C9290","#164B68","#D88418"])


def read(path):
    with path.open(newline="") as f: rows=list(csv.DictReader(f))
    return {k:np.array([float(r[k]) for r in rows]) for k in rows[0]}


def label(ax,s,x=-.10,y=1.04):
    ax.text(x,y,s,transform=ax.transAxes,fontsize=9.8,fontweight="bold",
            ha="left",va="bottom",color=INK,clip_on=False)


def bare(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)


def tank(ax):
    """Ideal Bloch-wave readout; finite-beam details are reserved for Fig. 5."""
    bare(ax); label(ax,"a",x=-.04,y=1.01)
    ax.set_title("Bloch-wave scattering",pad=4)
    ax.set_xlim(0,1); ax.set_ylim(0,1)
    ax.add_patch(Rectangle((.04,.08),.92,.79,fc=WATER,ec="#B8D4D7",lw=.75))

    # Periodic rigid metagrating.  The repeated unit cell is deliberately
    # schematic here; Fig. 2 carries the dimensional geometry.
    y=.24; left=.16; right=.90
    ax.plot([left,right],[y,y],color=INK,lw=1.0)
    period=(right-left)/7
    for m in range(7):
        x=left+m*period
        ax.plot([x+.18*period,x+.18*period,x+.50*period,x+.50*period],
                [y,y-.105,y-.105,y],color=INK,lw=.65)
        ax.plot([x+.64*period,x+.64*period,x+.76*period,x+.76*period],
                [y,y-.050,y-.050,y],color=INK,lw=.65)
    ax.annotate("",xy=(left+period,y-.14),xytext=(left,y-.14),
                arrowprops=dict(arrowstyle="<->",color=GREY,lw=.65))
    ax.text(left+.5*period,y-.19,r"$a$",ha="center",color=GREY)
    ax.text(.08,.25,r"$\cdots$",ha="center",va="center",fontsize=10,color=GREY)
    ax.text(.95,.25,r"$\cdots$",ha="center",va="center",fontsize=10,color=GREY)

    # Incident Bloch wave and the two reflected orders used in the response
    # maps.  The arrows terminate at the same aperture point to make the
    # channel projection visually unambiguous.
    hit=(.56,.34)
    ax.add_patch(FancyArrowPatch((.28,.73),hit,arrowstyle="-|>",mutation_scale=8,
                                 lw=1.05,color=GOLD))
    ax.text(.24,.77,r"$p_{\mathrm{in}}$",color=GOLD,ha="center")
    ax.add_patch(FancyArrowPatch(hit,(.72,.70),arrowstyle="-|>",mutation_scale=8,
                                 lw=1.05,color=BLUE))
    ax.add_patch(FancyArrowPatch(hit,(.82,.34),arrowstyle="-|>",mutation_scale=8,
                                 lw=1.05,color=GOLD))
    ax.text(.75,.74,r"$r_0$",color=BLUE,ha="center")
    ax.text(.86,.30,r"$r_{-1}$",color=GOLD,ha="center")

    # A phase-resolved line scan is the observable that produces the Fourier
    # amplitudes; it is not an additional finite-sample schematic.
    ax.plot([.21,.79],[.53,.53],color=TEAL,lw=.8,ls=(0,(3,2)))
    for xx in np.linspace(.25,.75,6):
        ax.add_patch(Circle((xx,.53),.007,fc="white",ec=TEAL,lw=.65))
    ax.text(.50,.58,r"$p(x)\;\longrightarrow\;r_n(k_x,f)$",
            ha="center",color=TEAL)


def response_data():
    r=read(DATA/"adaptive_scattering_response.csv")
    k=np.unique(r["kappa"]); u=np.unique(r["normalized_detuning"])
    nk,nu=len(k),len(u)
    def mat(name): return r[name].reshape((nk,nu),order="F")
    th=mat("theta_deg")[:,0]; fp=mat("pole_frequency_hz")[:,0]
    gam=mat("linewidth_hz")[:,0]
    freq=mat("frequency_hz")
    a0=mat("A0_real")+1j*mat("A0_imag")
    am=mat("Am1_real")+1j*mat("Am1_imag")
    bg0=.5*(a0[:,0]+a0[:,-1]); bgm=.5*(am[:,0]+am[:,-1])
    return th,fp,gam,freq,a0-bg0[:,None],am-bgm[:,None]


def channel_maps(fig,spec,th,fp,gam,freq,d0,dm):
    outer=fig.add_subplot(spec); bare(outer); label(outer,"b")
    outer.set_title("Simultaneous linewidth collapse",pad=4)
    sg=spec.subgridspec(2,1,hspace=.10)
    a0=fig.add_subplot(sg[0]); am=fig.add_subplot(sg[1],sharex=a0)
    y=np.linspace(-3.2,3.2,641)
    maps=[]
    global_scale=max(np.nanmax(abs(d0)),np.nanmax(abs(dm)))
    for data in (d0,dm):
        zz=np.full((len(th),len(y)),-4.0)
        for i in range(len(th)):
            good=np.isfinite(data[i].real)
            if good.sum()<3: continue
            x=freq[i,good]-fp[i]; val=abs(data[i,good])/global_scale
            order=np.argsort(x); x=x[order]; val=val[order]
            assert np.all(np.diff(x)>=0)
            z=np.interp(y,x,val,left=1e-4,right=1e-4)
            safe=np.clip(z,1e-4,None)
            zz[i]=np.log10(safe)
        maps.append(zz)
    assert np.all(np.diff(th)>0)
    thd=np.linspace(th.min(),th.max(),401); gamd=np.interp(thd,th,gam)
    dense=[]
    for zz in maps:
        dense.append(np.vstack([np.interp(thd,th,zz[:,j]) for j in range(len(y))]).T)
    for ax,zz,txt,col in [(a0,dense[0],r"$n=0$",TEAL),(am,dense[1],r"$n=-1$",GOLD)]:
        im=ax.pcolormesh(thd,y,zz.T,shading="auto",cmap=LINE,vmin=-4,vmax=0,rasterized=True)
        ax.plot(thd,.5*gamd,color="white",lw=.72); ax.plot(thd,-.5*gamd,color="white",lw=.72)
        ax.text(.015,.68,txt,transform=ax.transAxes,color=col,fontweight="semibold")
        ax.set_xlim(th.min(),th.max()); ax.set_ylim(y.min(),y.max())
        ax.set_ylabel(r"$f-\mathrm{Re}\,f_p$ (Hz)")
    a0.tick_params(labelbottom=False)
    am.set_xlabel(r"Angle $\theta$ (deg)")
    a0.scatter([th[0]],[0],s=21,c=RED,ec="white",lw=.45,zorder=5)
    a0.annotate("Rayleigh BIC",xy=(th[0],0),xytext=(th[0]+.10,2.45),
                color=RED,ha="left",va="top",
                arrowprops=dict(arrowstyle="-",color=RED,lw=.65))
    cb=fig.colorbar(im,ax=[a0,am],fraction=.027,pad=.018)
    cb.set_label(r"$\log_{10}$ resonant contrast"); cb.set_ticks([-4,-2,0]); cb.ax.tick_params(labelsize=7.3)


def linecuts(ax,th,fp,gam,freq,d0):
    label(ax,"c"); ax.set_title("Frequency cuts",pad=4)
    targets=[.12,.45,1.05]; cols=[RED,TEAL,BLUE]
    for dt,col in zip(targets,cols):
        i=int(np.argmin(abs((th-th[0])-dt)))
        x=freq[i]-fp[i]; y=abs(d0[i]); good=np.isfinite(y)
        x=x[good]; y=y[good]/max(y[good])
        order=np.argsort(x); x=x[order]; y=y[order]
        ax.plot(x,y,color=col,lw=1.35,label=rf"$\Delta\theta={th[i]-th[0]:.2f}^\circ$")
        # Open markers are exact theory subsamples reserved for experiment.
        ids=np.linspace(0,len(x)-1,13).round().astype(int)
        ax.scatter(x[ids],y[ids],s=13,fc="white",ec=col,lw=.55,zorder=4)
    ax.set_xlim(-3.2,3.2); ax.set_ylim(-.02,1.08)
    ax.set_xlabel(r"$f-\mathrm{Re}\,f_p$ (Hz)"); ax.set_ylabel("Normalized contrast")
    ax.legend(loc="upper right",fontsize=7.3,handlelength=1.4)


def extracted(fig,spec,summary):
    outer=fig.add_subplot(spec); bare(outer); label(outer,"d")
    outer.set_title("Extracted observables",pad=4)
    sg=spec.subgridspec(2,1,hspace=.12)
    ag=fig.add_subplot(sg[0]); aw=fig.add_subplot(sg[1],sharex=ag)
    dt=summary["theta_deg"]-summary["theta_deg"][0]
    g=np.maximum(summary["linewidth_hz"],1e-7)
    ag.semilogy(dt,g,color=PURPLE,lw=1.4)
    ids=np.linspace(2,len(dt)-2,12).round().astype(int)
    ag.scatter(dt[ids],g[ids],s=13,fc="white",ec=PURPLE,lw=.55,zorder=4)
    ag.set_ylabel(r"$\Gamma$ (Hz)"); ag.tick_params(labelbottom=False)
    w0=np.maximum(summary["spectral_weight_A0"],1e-10)
    wm=np.maximum(summary["spectral_weight_Am1"],1e-10)
    aw.semilogy(dt,w0,color=TEAL,lw=1.35,label=r"$n=0$")
    aw.semilogy(dt,wm,color=GOLD,lw=1.35,label=r"$n=-1$")
    aw.scatter(dt[ids],w0[ids],s=12,fc="white",ec=TEAL,lw=.5)
    aw.scatter(dt[ids],wm[ids],s=12,fc="white",ec=GOLD,lw=.5)
    aw.set_ylabel(r"$\int|\Delta A_n|^2df$"); aw.set_xlabel(r"$\theta-\theta_{\rm BIC}$ (deg)")
    aw.legend(loc="lower right",ncol=2,fontsize=7.3,handlelength=1.3)
    for a in (ag,aw): a.set_xlim(0,dt.max()); a.tick_params(labelsize=7.3,pad=1.6)


def field(fig,spec):
    outer=fig.add_subplot(spec); bare(outer); label(outer,"e")
    outer.set_title("Near-field reconstruction",pad=4)
    sg=spec.subgridspec(1,2,width_ratios=[1.55,.78],wspace=.27)
    af=fig.add_subplot(sg[0]); ak=fig.add_subplot(sg[1])
    x=np.loadtxt(DATA/"field_x_over_a.csv",delimiter=",")
    y=np.loadtxt(DATA/"field_y_over_a.csv",delimiter=",")
    pr=np.loadtxt(DATA/"field_total_real.csv",delimiter=",")
    sr=np.loadtxt(DATA/"field_scattered_real.csv",delimiter=",")
    si=np.loadtxt(DATA/"field_scattered_imag.csv",delimiter=",")
    meta=read(DATA/"field_metadata.csv")
    far_mask=y>=.8
    field_scale=np.nanpercentile(abs(sr[far_mask]),99.0)
    af.imshow(sr/field_scale,origin="lower",aspect="auto",extent=[x.min(),x.max(),y.min(),y.max()],
              cmap=FIELD,vmin=-1,vmax=1,interpolation="bilinear",rasterized=True)
    # Overlay sparse scan locations without hiding the continuous field.
    ix=np.arange(20,len(x),55); iy=np.arange(10,len(y),42)
    xx,yy=np.meshgrid(x[ix],y[iy]); af.scatter(xx,yy,s=4,fc="none",ec=INK,lw=.25,alpha=.45)
    af.set_xlabel(r"$x/a$"); af.set_ylabel(r"$y/a$"); af.set_xlim(-5,5); af.set_ylim(.12,3.2)
    p=sr+1j*si; jy=int(np.argmin(abs(y-1.5))); win=np.hanning(len(x))
    sp=np.fft.fftshift(np.fft.fft(p[jy]*win)); dx=x[1]-x[0]
    kxa=-2*np.pi*np.fft.fftshift(np.fft.fftfreq(len(x),d=dx))
    # k0*a from the final physical period retained by the model.
    omega=meta["frequency_hz"][0]*7.416665789e-3/1500
    q=kxa/(2*np.pi*omega); amp=abs(sp)/max(abs(sp)); mask=(q>-1.5)&(q<1.5)
    order=np.argsort(q[mask]); ak.plot(q[mask][order],amp[mask][order],color=BLUE,lw=1.2)
    theta=np.deg2rad(meta["theta_deg"][0]); q0=np.sin(theta); qm=q0-1/omega
    for qn,col,n in [(q0,TEAL,"0"),(qm,GOLD,"−1")]:
        ak.axvline(qn,color=col,lw=.9,ls=(0,(2,2))); ak.text(qn,.94,n,color=col,ha="center",va="top")
    ak.set_xlim(-1.5,1.5); ak.set_ylim(0,1.03); ak.set_xlabel(r"$k_x/k_0$"); ak.set_ylabel("Spatial spectrum")


def main():
    th,fp,gam,freq,d0,dm=response_data(); summary=read(DATA/"pole_and_radiation_summary.csv")
    assert np.all(np.diff(th)>=0)
    fig=plt.figure(figsize=(7.15,5.65),facecolor="white")
    gs=fig.add_gridspec(2,6,left=.06,right=.985,bottom=.09,top=.95,
                        height_ratios=[.92,1.08],wspace=.75,hspace=.40)
    tank(fig.add_subplot(gs[0,0:2]))
    channel_maps(fig,gs[0,2:6],th,fp,gam,freq,d0,dm)
    linecuts(fig.add_subplot(gs[1,0:2]),th,fp,gam,freq,d0)
    extracted(fig,gs[1,2:4],summary)
    field(fig,gs[1,4:6])
    OUT.mkdir(parents=True,exist_ok=True); stem=OUT/"fig4_scattering_observability_nature"
    fig.savefig(stem.with_suffix(".pdf"),bbox_inches="tight",facecolor="white")
    fig.savefig(stem.with_suffix(".svg"),bbox_inches="tight",facecolor="white")
    fig.savefig(stem.with_suffix(".png"),dpi=600,bbox_inches="tight",facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"),dpi=600,bbox_inches="tight",facecolor="white",
                pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig); print(stem)


if __name__=="__main__": main()
