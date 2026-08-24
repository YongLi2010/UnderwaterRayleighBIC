# Reproduction audit of the present Rayleigh-BIC manuscript

Date: 2026-08-24  
Environment: MATLAB R2026a Update 3 (`26.1.0.3276743`) on macOS  
Scope: read-only audit of `Ni2019_MATLAB`, the cached `.mat/.csv` results, and
the current `arxiv_theory_paper/main.tex` and `suppl.tex`.  No MATLAB source,
TeX source, bibliography, or figure file was edited.  The only files written
by this audit are this report and raw logs under
`advanced_analysis/reproduction_raw/`.

## Executive summary

The final manufacturable strict-root calculation is independently reproducible
at `(N,K)=(313,39)`.  The direct rerun gives

```text
kappa = 0.08384869655290064
Omega = 0.91615130344709939
theta = +/-5.251216456430976 deg
sigma_min/sigma_max = 1.2603643174303456e-15
r_raw = 1.3765667611238457e-13
```

The dimensional geometry and the homogeneous radiation phasor cancellation
also reproduce the manuscript.  The off-Rayleigh point and the reciprocal
two-port audit reproduce to the quoted precision.

Two qualifications are important before submission:

1. The pole cache used by the current figure generator is not the final
   `(313,39)` strict-root geometry.  It is the `(N,K)=(121,15)` intermediate
   root with `kappa=0.08394063684879589`, whose geometry differs from the final
   root by up to `2.01e-4` in normalized coordinates.  The pole cache itself
   is exactly reproducible, including `Q_max=2.0265381467721e12` and the local
   slope `-1.99608687`, but this should be identified as a finite-truncation
   continuation diagnostic rather than silently presented as a continuation
   of the final `(313,39)` root.

2. The router-convergence table in `suppl.tex` and the hard-coded values in
   `run_PRL_candidate_figure_gallery.m` are not reproduced at their stored
   frequencies for the final geometry.  Direct values at `(121,15)`,
   `(185,23)`, and `(281,35)` are `5.43e-4`, `5.04e-3`, and `0.1469`, rather
   than `0.30737`, `0.30662`, and `0.30620`.  The final `(313,39)` point gives
   `0.306061`, close to the stated `0.30613`.  This table must be regenerated
   from a documented frequency-retuning procedure or removed.

The default strict-operator threshold classifier also sees the target
`k_y=1.1921e-7` as a tiny positive open value because the default absolute
tolerance is `5.7563e-8`.  The target amplitude is nevertheless explicitly
removed, and the residual is unchanged.  With `KyTolerance=1e-6`, the same
calculation is classified as a threshold.  This is a numerical classification
detail, not a change in the normalized Rayleigh condition
`Omega=|kappa-1|`.

## Model and normalization traced from source

The driven solver is `Ni2019_MATLAB/ni2019_modal_solver.m`:

- `lambda`, `a`, groove widths/depths/gap, and all coordinates use one common
  length unit;
- the external orders are `n=-(N-1)/2,...,(N-1)/2`;
- each groove retains cosine modes `q=0,...,K-1` with a rigid bottom;
- time dependence is `exp(+j omega t)`;
- the reflected field uses `exp(-j k_y y)` and the outgoing/decaying branch;
- driven efficiency is
  `eta_n=Re(k_y,n)/k_y,inc * abs(A_n/A_inc)^2`, with `A_inc=1`.

The homogeneous pole-free operator is
`Ni2019_MATLAB/ni2019_full_eigen_operator.m`.  The strict operator is
`Ni2019_MATLAB/ni2019_strict_rayleigh_operator.m`; it removes the selected
threshold amplitude and every finite-flux open amplitude before scaling and
SVD.  The reduced physical mode is reconstructed from the scaled right
singular vector, and the removed amplitudes are inserted as exact zeros.

For the final dimensional design, the cache
`Ni2019_MATLAB/results/StrictRayleighBIC_200kHz_min1mm.mat` stores

```text
xFinal = [kappa,d1/a,d2/a,w1/a,w2/a,g/a]
       = [0.08384869655290064,
          0.66094150824509523,
          0.15726373197723323,
          0.57701824931554235,
          0.15325679426861946,
          0.13486247820791922]
```

At `f=200 kHz`, `c=1500 m/s`,
`a=Omega*c/f=6.87113477585325 mm` and `lambda=c/f=7.5 mm`.  The resulting
dimensions are

```text
d1 = 4.54141818210777 mm       d2 = 1.08058029776923 mm
w1 = 3.96477015917398 mm       w2 = 1.05304808873490 mm
g  = 0.926658263972184 mm
```

Both groove widths satisfy the stated `1 mm` constraint.  The inter-groove
gap is about `0.927 mm`; the manuscript does not claim that the gap is at
least `1 mm`.

## Claim-by-claim verification table

`Rerun` means the quantity was recomputed in this audit by calling the named
routine in memory.  `Cache` means it was read from an existing `.mat/.csv`.
The cache was never overwritten.

| Manuscript quantity | Source routine/cache | Independent result | Verdict |
|---|---|---|---|
| Final BIC geometry, `kappa`, `Omega`, `theta` | `run_Ni2019_strict_rayleigh_bic_200k_1mm.m`; `StrictRayleighBIC_200kHz_min1mm.mat/.csv` | Exact values above; `theta=5.251216456430976 deg` | Reproduced |
| Dimensional period and dimensions | same cache, direct conversion | `a=6.87113477585325 mm`; all values above | Reproduced |
| Strict reduced SVD at `(313,39)` | direct `ni2019_strict_rayleigh_operator` | `sigma_rel=1.2603643174303456e-15`; `r_raw=1.3765667611238457e-13` | Reproduced |
| Full pole-free operator at exact normalized endpoint | direct `ni2019_full_eigen_operator_complex` with `kappa=xFinal(1)`, `Omega=1-kappa`, and target `k_y=0`, evaluated on the strict mode | full-operator `sigma_rel=7.715113792404351e-16`; strict-mode raw residual `1.3765667768908454e-13` | Reproduced at the final `(313,39)` truncation |
| Negative-angle partner | direct strict operator with `theta=-5.251216456430976 deg`, target `n=+1` | `sigma_rel=1.201286449663739e-15`; `r_raw=1.105241300315003e-13`; removed orders `0,+1` | Reproduced, with expected finite-truncation asymmetry |
| Homogeneous selected amplitudes | direct strict operator | positive-angle `A_-1=0`, `A_0=0` exactly after constrained reconstruction; `A_+1=+0.23701341906568918 i` in the chosen mode normalization | Reproduced |
| Radiation phasor cancellation | direct velocity-block decomposition of final strict mode | `A_0` groove contributions sum to `7.92e-17`; `A_-1` contributions sum to `5.00e-17` (relative `8.36e-16` and `5.18e-16`) | Reproduced |
| Evanescent exterior localization | direct strict-mode energy proxy | `E_ext(H=1a)=0.0176901497887`; saturates at `0.0177019452335` by `H=5a` and remains constant to shown precision at `20a` | Reproduced |
| Rayleigh classification | direct `ni2019_full_eigen_operator` | normalized `Omega=1-kappa` exactly, but floating-point `k_y,-1=1.1920928955e-7`; default strict tolerance `5.7563484090e-8`, so `target_is_rayleigh=false` and `threshold_orders=[]` | Numerical caveat; use a documented larger threshold tolerance or direct normalized `kappa/Omega` construction |
| Cached pole endpoint | `StrictRayleighBIC_200kHz_min1mm_poles.mat`; direct rerun of `ni2019_track_leaky_pole_to_rayleigh` with cached cfg | cache cfg is intermediate root `kappa0=0.08394063684879589`, `Omega0=0.91605936315120406`, `N/K=121/15`; rerun gives `same_branch=1`, `Q_max=2.026538146772098e12`, `q_end=-4.64217e-11-4.46005776e-3 i`, `|Omega_end-Omega0|=8.57384e-7` | Cache and rerun agree |
| Local pole slope | same pole cache/rerun | over `3e-5 <= |Delta kappa| <= 1e-3`, slope `-1.99608687`; all 19 finite points slope `-1.99456136` | Reproduced as a finite-model diagnostic |
| Pole relation to final strict root | direct exact-final geometry diagnostic with `N/K=121/15` | endpoint `sigma_rel=4.095e-4`, `q_end=-1.45381e-8-1.30981e-2 i`, `|Omega_end-Omega0|=8.3635e-5` at `|Delta kappa|=1e-5` | Not the same numerical continuation; clarify or recompute |
| Practical router, fixed final geometry | direct `ni2019_modal_solver`, `theta=+6 deg`, `f=200216.874212 Hz`, `N/K=313/39` | `eta_-1=0.306061126315771`; `eta_0=0.693938873684255`; `theta_out,-1=-80.336964032 deg`; energy error `2.665e-14` | Numerically reproduced; manuscript rounding differs by `6.69e-5` in `eta_-1` |
| Reverse router | same direct solver, `theta=-6 deg` | `eta_+1=0.306061126316284`; `eta_0=0.693938873684022`; `theta_out,+1=+80.336964032 deg` | Reproduced |
| Router convergence table | hard-coded `routeConvergence...` in `run_PRL_candidate_figure_gallery.m:62-65` and `suppl.tex:242-245` | Direct fixed-final-geometry values at stored frequencies: `(121,15)=0.0005427505`, `(185,23)=0.005038793`, `(281,35)=0.1468842`, `(313,39)=0.3060611` | **Not reproduced; high-priority correction** |
| Ideal hot point | cache `PRL_article_figure_data_v5.mat`; direct final-geometry solver | `theta=5.270442546972 deg`, `f=200005.779468 Hz`, `(313,39)` gives `eta_-1=0.8568063583`, `eta_0=0.1431936432`, `cond=1.6209834703e10`, energy error `1.45e-9` | Reproduced; explicitly ill-conditioned and not robust |
| Off-Rayleigh state | direct `ni2019_modal_solver` | `theta=32.6328 deg`, `f=202430 Hz`, `Omega=0.927282541783982`, `n=-1` threshold `141824.999447620 Hz`; at `(401,50)`, `eta_-1=0.999618503882144`, `eta_0=0.000381496117855`, closure error `4.441e-16` | Reproduced |
| Off-Rayleigh convergence | same direct solver | `(81,11): .9779112`; `(121,15): .9929755`; `(161,20): .9986358`; `(201,25): .9996289`; `(313,39): .9999151`; `(401,50): .9996185` | Reproduced; nonmonotonic last step is visible and should remain disclosed |
| Single-cavity ablation | `SingleGroove_strict_Rayleigh_search.mat/.csv`; direct fixed-final-geometry strict operator | cached reoptimized residual at `(221,44)` `sigma_rel=.1245882`, raw residual `1.48399`; direct cross at `(241,48)` `sigma_rel=.1236552`, raw residual `1.49098` | “About 0.124” is reproduced as a singular-value ratio; not a universal one-cavity no-go |
| Reciprocity/unitarity | direct power-normalized matrix reconstruction | at route `(313,39)`: reciprocity Frobenius error `1.97401310518e-12`; unitarity errors `3.91671503963e-13` and `3.72208418564e-12` | Reproduced |

## Direct driven route data at fixed frequency

The following table is useful because it exposes the truncation sensitivity that
the current hard-coded “router convergence” table hides.  All rows use the
same final geometry, `theta=+6 deg`, and `f=200216.874212 Hz`; only `(N,K)`
changes.

| `(N,K)` | `eta_-1` | `eta_0` | condition number | energy error |
|---:|---:|---:|---:|---:|
| `(81,11)` | `9.1843e-5` | `0.9999082` | `1.954e4` | `2.2e-16` |
| `(121,15)` | `4.0379e-4` | `0.9995962` | `4.680e4` | `1.4e-15` |
| `(201,25)` | `6.5929e-3` | `0.9934071` | `2.610e5` | `8.9e-16` |
| `(313,39)` | `0.3060611` | `0.6939389` | `2.593e6` | `2.7e-14` |
| `(401,50)` | `0.0571528` | `0.9428472` | `1.408e6` | `1.2e-14` |

This is not a contradiction of a finite-truncation route point, but it means
that a publication-quality claim must specify whether frequency and geometry
are retuned jointly at every `(N,K)`, and must show the resulting retuned
frequencies.  The four frequencies hard-coded in the present gallery do not
produce the four stored efficiencies when passed to the current solver.

## Additional requested physics not present as reproducible results

The current code/manuscript set does not contain a completed, independently
rerunnable version of the following requested analyses:

- a `c=1400...1600 m/s` environmental sweep with a table of BIC frequency
  shifts and radiation residuals;
- dimensional scale-factor and `f_BIC(c)` plots with normalized parameters
  held fixed;
- random independent perturbations of all five dimensional variables
  `d1,d2,w1,w2,g` at `+/-1%`, `+/-2%`, and `+/-5%`, with distributions of
  homogeneous residual and threshold amplitude error;
- a three-way homogeneous-field comparison of strict BIC, near-BIC, and an
  ordinary resonance using one consistent eigenfield normalization.

There is a cached fixed-root 2D map varying only `d2` and `w2` in
`PRL_article_figure_data_v5.mat`, and the Supplement currently describes the
sound-speed scaling as a guide (`1500 -> 1480 m/s`, predicted retuning near
`197.33 kHz`).  Those are not substitutes for the requested 1400--1600 m/s
sweep or five-parameter random tolerance study.

## Validation commands and raw logs

All solver calls below were in-memory calls and did not invoke manuscript or
figure-writing runners:

```text
01_strict_operator.txt                 final/negative strict SVD and zeros
02_cache_values.txt                    v5 cache values used by figures
03_driven_router_corrected.txt         driven route, fixed frequency, N/K sweep
04_off_rayleigh.txt                    off-Rayleigh direct rerun
05_pole_cache.txt                      cached pole fields and trajectory
05b_pole_slope_cfg.txt                 cached pole cfg metadata
06_pole_rerun_final_geometry.txt       exact-final pole diagnostic
07_pole_rerun_cached_cfg.txt           exact cached-config pole rerun
08_router_retuning.txt                 hard-coded table versus direct rerun
09_single_cavity.txt                   cached and direct single-cavity checks
10_refine_final_root.txt               independent in-memory final-root polish
11_threshold_tolerance.txt             target-threshold classifier tolerance
12_hot_point_corrected.txt             hot quasi-BIC and practical route checks
13_reciprocity_direct.txt              direct S-matrix reciprocity/unitarity
14_phasors.txt                         two-groove radiation phasors
15_exterior_energy.txt                 evanescent energy saturation
16_full_operator_final_endpoint.txt    final full-operator endpoint check
```

The successful pole slope was also independently evaluated from the cached
numeric arrays with SciPy; it gives `-1.99608687` on the manuscript's stated
finite interval and `-1.99456136` over all 19 finite cached points.  The raw
logs from exploratory failed command syntaxes are retained only as audit
history and are not used as validation evidence.

## Recommended decisions for the primary manuscript thread

Before submission, the primary thread should decide whether to:

1. recompute a pole continuation at the final `(313,39)` root (or explicitly
   relabel the current pole plot as an intermediate-truncation diagnostic);
2. regenerate the router-retuning table from an actual optimizer and record
   every frequency, `(N,K)`, condition number, and energy error;
3. state the threshold classifier tolerance or evaluate the strict operator
   using an exact normalized `(kappa,Omega)` input;
4. add the missing sound-speed, scale, random five-parameter tolerance, and
   homogeneous-field studies if they are to remain part of the submission
   package.

No claim of experimental validation is supported by this audit.
