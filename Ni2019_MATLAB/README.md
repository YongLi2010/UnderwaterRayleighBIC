# Ni2019 acoustic meta-grating reproduction

This folder provides a general modal-matching program and reproduces every
distinct geometry printed in Figs. 2-6 of:

H. Ni, X. Fang, Z. Hou, Y. Li, and B. Assouar, "High-efficiency anomalous
splitter by acoustic meta-grating," *Physical Review B* **100**, 104104
(2019), DOI: 10.1103/PhysRevB.100.104104.

## Run

Open MATLAB in this folder and execute the comprehensive comparison:

```matlab
run_Ni2019_all_results
```

This writes CSV, MAT, and PNG results under `results/`. For the original
four-case quick run, execute:

```matlab
run_Ni2019_reproduction
```

The script prints every propagating diffraction order, its angle, and its
normal-power efficiency

```text
eta_n = Re(ky_n)/ky_inc * abs(A_n/A_inc)^2
```

It also produces an efficiency bar chart and a truncation-convergence plot.
The returned structures `Ni2019_cases` and `Ni2019_results` are placed in the
base workspace for further analysis.

## Method and conventions

- External field: Floquet diffraction expansion.
- Groove field: cosine waveguide modes satisfying the rigid groove bottom.
- Interface: pressure and normal-velocity continuity projected in a Galerkin
  modal-matching formulation equivalent to Appendix A of the paper.
- Time convention: `exp(+j*omega*t)`.
- Default truncation: `N=101` external modes and `K=10` modes per groove,
  matching the paper.
- The groove group is centered within a period. A common lateral translation
  only changes diffraction phases and does not change efficiencies.
- The calculation is lossless. Thermoviscous loss, finite sample width, source
  aperture, and experimental boundaries are outside this first reproduction.

For the numerical pressure amplitudes printed in the article, the complete run
gives a mean relative difference of `0.1925%` and a maximum of `0.7505%`.
Lossless energy closure is typically at the `1e-15` level.

## Important 5:3 convention check

For the geometry labeled `I_-1:I_+1 = sqrt(5):3` in the paper
(`t/a=[0.159,0.230]`, `d/lambda=[0.229,0.170]`, `dx/a=0.176`),
the solver reproduces the paper's reported pressure amplitudes: their ratio is
approximately `5:3`. Because the two beams leave at symmetric angles, however,
the normal-power ratio is the squared amplitude ratio, approximately `25:9`.
The computed efficiencies are about `0.734` and `0.266`, revealing an internal
inconsistency between the paper's power-ratio label and its quoted amplitudes.

Two more printed-parameter diagnostics are retained rather than hidden:

- The Fig. 4(a) geometry printed for a 0:1 target gives efficiencies
  `0.0987, 0.4043, 0.4970` in orders `-1,0,+1`. The optional optimizer finds a
  separate valid 0:1 design but does not claim that it is a correction.
- The Fig. 5 period and grating equation give `-29.594 deg` for order -1, while
  the text says `-28.4 deg`. Its printed amplitudes and equal-power split are
  nevertheless reproduced accurately.

## Files

- `run_Ni2019_reproduction.m`: published cases, report, plots, convergence.
- `ni2019_modal_solver.m`: reusable modal-matching solver.
- `ni2019_case_library.m`: all six geometries printed in the paper.
- `run_Ni2019_all_results.m`: complete calculation, comparison table, and figures.
- `ni2019_reconstruct_field.m`: pressure, velocity, and intensity reconstruction.
- `ni2019_optimize_geometry.m`: optional bounded geometry optimizer.
- `run_Ni2019_reoptimize_F4a.m`: independently re-designs the inconsistent
  printed Fig. 4(a) case; it does not silently alter the paper's parameters.
- `ni2019_full_eigen_operator.m`: pole-free homogeneous operator retaining
  both outgoing Floquet and groove-mode unknowns.
- `ni2019_find_rayleigh_bic.m`: two-depth search at `a=lambda0` using the
  smallest singular value and open-channel radiation as joint criteria.
- `run_Ni2019_rayleigh_bic_depth_study.m`: depth maps, scattering threshold
  diagnostics, and depth/frequency Argand trajectories.
- `ni2019_full_eigen_operator_complex.m`: complex-frequency operator with an
  explicitly selected Rayleigh channel and Riemann sheet.
- `ni2019_continue_complex_pole.m`: complex-q pole continuation versus Bloch
  wavenumber.
- `run_Ni2019_complex_pole_Q_scaling.m`: paper-colored pole, Q-scaling, and
  Riemann-sheet plots.
- `run_Ni2019_reciprocal_direction_pair.m`: compares +k and -k poles and
  demonstrates their time-reversal pairing.
- `ni2019_optimize_single_rayleigh_bic.m`: constrained multistart search for
  an off-Gamma, single-Rayleigh threshold state.
- `ni2019_continue_offgamma_pole.m`: physical/improper-sheet pole continuation
  about an arbitrary off-Gamma Rayleigh endpoint.
- `run_Ni2019_single_rayleigh_bic_validation.m`: historical first-stage
  threshold-pole diagnostic; the stricter test below supersedes its BIC label.
- `run_Ni2019_strict_rayleigh_bic_verification.m`: applies three stricter
  checks: independent leaky-pole continuation, a separately evaluated
  aperture-radiation numerator, and square-integrability/energy diagnostics.
- `ni2019_strict_rayleigh_operator.m`: removes the target grazing amplitude
  and every finite-flux open amplitude before solving the homogeneous system.
- `ni2019_optimize_strict_rayleigh_bic.m`: particle-swarm/multistart discovery
  of strict zero-radiation basins over several modal truncations.
- `ni2019_refine_strict_rayleigh_bic.m`: variable-projection least-squares
  polishing of a strict compatibility zero at one finite truncation.
- `run_Ni2019_strict_rayleigh_bic_optimization.m`: practical increasing-order
  root sequence, cross-truncation check, complex-pole track, and paper-colored
  summary figure.
- `run_Ni2019_single_rayleigh_scattering_characterization.m`: driven
  angle-frequency maps of `A_-1,A_0,A_+1`, fixed-angle/fixed-frequency
  amplitude-phase-efficiency cuts, and representative pressure/intensity
  fields for experimental characterization near the strict BIC.
- `run_Ni2019_single_groove_strict_search.m`: wide-aperture single-groove
  ablation/search retaining up to 44 transverse groove modes; performs a
  global survey, increasing-order reoptimization, and fixed-geometry
  cross-truncation comparison against the strict two-groove root.
- `run_Ni2019_strict_rayleigh_bic_200k_1mm.m`: constrained two-groove
  reoptimization at 200 kHz in nominal 1500 m/s water, enforcing a guaranteed
  physical groove-width minimum of 1 mm and reporting the resulting geometry
  directly in millimeters.
- `ni2019_optimize_three_groove_double_rayleigh_bic.m`: global/multistart
  full or mirror-symmetric three-groove search at the Gamma double threshold.
- `ni2019_refine_three_groove_double_rayleigh_bic.m`: eight-variable
  variable-projection polishing with `A_-1=A_0=A_+1=0` imposed directly.
- `ni2019_three_groove_parity_diagnostic.m`: reflection-parity and internal
  groove-mode decomposition for the symmetric three-groove eigenmode.
- `run_Ni2019_three_groove_double_rayleigh_bic_optimization.m`: converged
  odd-mode root sequence, evanescent-energy test, field, pole/Q trend, and
  200-kHz dimensional design.
- `ni2019_two_groove_odd_strict_operator.m`: independent mirror-odd strict
  operator for two identical centered grooves at the Gamma double threshold.
- `ni2019_two_groove_even_strict_operator.m`: complementary mirror-even
  strict operator, including both the independent `n=0` and `n=+/-1`
  compatibility equations.
- `ni2019_two_groove_parity_diagnostic.m`: two-groove parity, threshold-source,
  and transverse-mode decomposition.
- `ni2019_optimize_two_symmetric_groove_double_rayleigh_bic.m`: global and
  multistart search in the two-groove mirror-odd sector.
- `ni2019_refine_two_symmetric_groove_double_rayleigh_bic.m`: odd-sector
  variable-projection root polishing at one finite truncation.
- `run_Ni2019_two_symmetric_groove_highmode_feasibility.m`: removes the wide
  center groove, tests both parities through `K=31`, compares the complete
  residual with the singular-value diagnostic, and plots the local `q=1`
  cutoff mechanism.

## Rayleigh-BIC depth study

### Strict two-groove optimization

For the current off-Gamma two-groove result, run:

```matlab
run_Ni2019_strict_rayleigh_bic_optimization
```

The strict operator imposes `A_-1=A_0=0` before the SVD; `A_+1` remains in
the model because it is evanescent at this single-Rayleigh point.  Increasing
finite truncations `(N,K)=(121,15),(185,23),(281,35),(313,39)` give polished
compatibility residuals between roughly `1e-13` and machine precision while
the design parameters stabilize at

```text
kappa       = 0.082502670831
Omega       = 0.917497329169
theta       = 5.159086970 deg
d/a         = [0.656894843544, 0.221773792742]
w/a         = [0.575896993065, 0.037613144177]
gap/a       = 0.193244931379
```

The maximum parameter change from `K=35` to `K=39` is `1.53e-5`; the groove
participation is about `0.983`.  The independently seeded outgoing-sheet pole
track reaches `Q>1e13` before double-precision noise dominates.  Over the
stable finite-model interval the observed fit is near `Q~|Delta kappa|^-2`.
This is reported as a numerical trend, not as a proof of the cubic law; a
separate reduced asymptotic model is still required to identify the ultimate
continuum exponent.

The default run starts from the documented basin seed so that it finishes
quickly.  Set the environment variable `NI2019_GLOBAL_SEARCH=1` before
starting MATLAB to repeat the particle-swarm discovery stage.  Set
`NI2019_SKIP_POLES=1` to run only the strict eigenvalue optimization.

The practical interpretation is a robust very-high-Q Rayleigh-BIC design.
Machine-infinite Q is not used as an engineering acceptance criterion; the
finite-root zero, parameter convergence, open-channel cancellation, and pole
trend are stored separately so their numerical roles are not conflated.

### Earlier candidates and diagnostics

Run `run_Ni2019_single_rayleigh_bic_search` followed by
`run_Ni2019_single_rayleigh_bic_validation` for the off-Gamma,
single-Rayleigh calculation. This validation is completed before any
unidirectional optimization is introduced.

The stricter three-part check currently classifies that optimized point as a
Rayleigh-threshold pole candidate rather than a verified BIC: the nontrivial
pole branch is confirmed, but the aperture-radiation zero is not stable under
groove-mode truncation and the nonzero grazing harmonic makes the exterior
energy grow linearly with integration height. Run
`run_Ni2019_strict_rayleigh_bic_verification` to reproduce this verdict.

### Three-groove Gamma double-Rayleigh BIC

The separate three-groove calculation fixes `kappa=0`, `Omega=1`, and removes
all three non-decaying exterior amplitudes `A_-1`, `A_0`, and `A_+1` before
solving the homogeneous system:

```matlab
run_Ni2019_three_groove_double_rayleigh_bic_optimization
```

The converged mirror-symmetric odd mode has

```text
d/a = [0.192681441083, 0.676102882473, 0.192681441083]
w/a = [0.091719249004, 0.529780736125, 0.091719249004]
g/a = [0.073658809190, 0.073658809190]
```

At `(N,K)=(313,39)`, its strict singular-value ratio is `6.7e-16`, the raw
compatibility residual is `9.1e-14`, and the maximum parameter change from
`K=31` to `K=39` is `7.3e-6`.  The mode is odd to a coefficient residual of
about `5.4e-12`: the two outer-groove `q=0` modes are opposite in phase and
the wide central groove is dominated by its `q=1` transverse mode.  This
explains why a three-groove fundamental-mode-only model misses the solution.

All propagating and threshold power is zero in the constrained eigenmode;
the remaining exterior harmonics are evanescent and their integrated energy
saturates with height.  Outgoing-sheet continuation over the numerically
stable interval gives `Q>2.9e12` and a local fit near
`Q~|kappa|^-2.04`.  The exponent is reported as a finite-model trend rather
than an imposed asymptotic law.

At 200 kHz in water with `c=1500 m/s`, the Gamma condition gives `a=lambda`
and therefore `a=7.5 mm`.  The generated CSV contains all normalized and
millimetre dimensions.  Set `NI2019_GAMMA3_GLOBAL_SEARCH=1` to repeat the
global discovery stage or `NI2019_GAMMA3_SKIP_POLES=1` to skip pole tracking.

### Three-groove search with 1 mm lateral manufacturing constraints

To require every groove width and every edge-to-edge gap (including the gap
across the periodic boundary) to be at least 1 mm at 200 kHz, run:

```matlab
run_Ni2019_three_groove_double_rayleigh_bic_min1mm
```

The independent mirror-parity operator
`ni2019_three_groove_parity_strict_operator` removes
`A_-1=A_0=A_+1=0` exactly and checks both odd and even sectors through
`(N,K)=(313,39)`.  The six lateral features share the 7.5 mm period, so the
center groove obeys

```text
w_center <= 7.5 mm - 5*(1 mm) = 2.5 mm.
```

The propagating center-groove `q=1` mode instead requires
`w_center >= a/2 = 3.75 mm`.  Consequently the original wide-center-groove
mechanism is geometrically closed: `beta_1^2*a^2 <= -5*pi^2` and
`gamma_1*a >= sqrt(5)*pi`.

Constrained symmetric odd/even and full asymmetric searches did not find a
converged strict root.  The best candidates move one or more widths to the
1 mm lower boundary; at the highest verification truncation their normalized
singular-value ratios remain of order `1e-2` and their complete unscaled
homogeneous residuals remain of order `1e-1`.  A literal `>1 mm` probe using
1.01 mm widths also fails.  Fresh odd/even optimization restarts with the
1.01 mm bound and depths extended through `3a` are saved separately;
neither hits the depth upper bound, but both remain at order `1e-2` in the
high-order singular-value test.  This is a numerical
no-root result for the searched rectangular three-groove design, not a
general mathematical no-BIC theorem.

### Three-groove strict root with all lateral features above 0.7 mm

Reducing the manufacturing limit to a strict `>0.70 mm` restores the wide
center-groove mechanism.  The reproducible study uses a practical numerical
minimum of 0.71 mm:

```matlab
run_Ni2019_three_groove_double_rayleigh_bic_min0p7mm
```

At 200 kHz (`a=7.5 mm`), the converged `K=95` finite-truncation geometry is

```text
depths       = [1.450216217, 5.341324975, 1.450216217] mm
widths       = [0.710003724, 3.949984714, 0.710003724] mm
internal gap = [0.710000000, 0.710000000] mm
boundary gap =  0.710007837 mm
```

The center width exceeds the local `q=1` cutoff `a/2=3.75 mm`, giving
`beta_1*a=1.973912`.  Sequentially repolishing only the outer and center
depths at `K=23,31,39,43,47,55,63,79,95` keeps the independent odd-parity
singular-value ratio near `3e-15` and the complete unscaled residual near
`2e-13`; the maximum parameter step continues to decrease.  Quadratic
extrapolation in `1/K` over the five highest roots gives continuum depths
`[1.450178784,5.341704350,1.450178784] mm`.  When that single extrapolated
geometry is held fixed, its singular-value ratio decreases monotonically
from `1.30e-4` at `K=39` to `2.44e-5` at `K=111`.  An independent generic
strict operator gives a singular-value ratio near `8e-16` at `K=95`, with all
three amplitudes `A_-1,A_0,A_+1` exactly removed and groove participation
about `0.989`.

Outgoing-sheet continuation from the finite `K=15` root gives 16 accepted
pole points, `Q>3.5e9` at the smallest accepted offset, and a local trend
`Q~|kappa|^-2.00`.  Set `NI2019_MIN0P7_REFINE=1` before running the script to
repeat the truncation-by-truncation depth polishing instead of using the
stored converged sequence.

### Two identical side grooves with higher transverse modes

The center-groove ablation is reproduced by

```matlab
run_Ni2019_two_symmetric_groove_highmode_feasibility
```

It retains the two outer grooves of the converged three-groove geometry,
replaces the central groove by rigid material, and tests the mirror-odd and
mirror-even sectors separately.  Both operators impose
`A_-1=A_0=A_+1=0` exactly.  The odd reduced matrix is one row taller than its
unknown vector; the even matrix is two rows taller because even parity leaves
both the `n=0` channel and one independent `n=+/-1` combination to cancel.

For the literal two-outer-groove clone,

```text
[d/a,w/a,g/a] = [0.192681441083,0.091719249004,0.677098354505].
```

At `(N,K)=(161,31)`, the odd and even singular-value ratios are approximately
`9.3e-2` and `1.4e-1`; their complete raw residuals remain about `0.33` and
`0.57`.  Broad searches of `d/a`, `w/a`, and `g/a`, including both parity
sectors and depths through `3a`, found only small-width/small-gap boundary
minima or modes dominated by the highest retained transverse orders.  Their
complete residual stays at order `1e-1` instead of converging to zero.  These
points are therefore radiation-dark or vanishing-aperture near-nulls, not
strict eigenmodes.

The geometric reason is visible from

```text
beta_1^2 a^2 = (2*pi)^2 - (pi/(w/a))^2.
```

Two equal non-overlapping grooves require `2w+g<a`, hence `w<a/2`; their
`q=1` mode is necessarily evanescent in depth at `Omega=1`.  For the retained
outer groove, `gamma_1*a=33.67` and `exp(-gamma_1*d)=1.52e-3`, so its apparent
high-order aperture content is a shallow boundary layer rather than an
independent cavity resonance.  By contrast, the converged three-groove center
has `w_c/a=0.52978>1/2`, a real `beta_1*a=2.07694`, and its exact mode is
dominated by the propagating `q=1` resonance.

This numerical result is deliberately stated as "no converged finite-geometry
root found," not as a formal no-go theorem.  A two-opening design can regain
the missing internal degree of freedom by using a widened/compound cavity
behind each narrow opening, a different cavity shape, or another independent
resonator; merely increasing `K` for two simple small rectangular grooves does
not create that resonance.

The study fixes the lateral geometry of the published -88 deg structure,
sets the exact first-order threshold `a=lambda0`, and varies only the two
groove depths:

```matlab
run_Ni2019_rayleigh_bic_depth_study
```

With `N=61` and `K=10`, the candidate is near
`d/lambda0=[0.48686,0.04077]`. The normalized smallest singular value is
about `4.7e-9`, the zeroth-order component of the null vector is below
`2e-11` in power fraction, and the two grazing amplitudes have equal
magnitudes near `0.321`. Increasing `N` and using `K>=6` moves the candidate
smoothly toward approximately `[0.487,0.040]`.

The eigenvalue calculation deliberately uses the full block operator. The
reduced scattering operator contains `tan(beta*d)` and therefore has poles at
quarter-wave groove resonances that can be mistaken for eigenvalue zeros.

Complex-pole continuation confirms `Im(Omega)->0`, but the Gamma-point design
contains simultaneous +1 and -1 Rayleigh branch points. Its fitted scaling is
approximately `Q~|kappa|^-0.65` (approaching a 2/3 fractional exponent), not
the cubic law of the single-Rayleigh off-Gamma model in 2607.07228v1. Unequal
depths break mirror symmetry but not time-reversal symmetry: the +k and -k
poles remain a reciprocal pair. A one-sided radiating state is therefore a
UGR/quasi-BIC unless the remaining channel is exactly grazing or closed.

For dimensional calculations, replace `lambda=1` by the desired wavelength and
keep all geometry in the same length unit. At the paper's 8 kHz air experiment,
`lambda=42.88e-3` m.
