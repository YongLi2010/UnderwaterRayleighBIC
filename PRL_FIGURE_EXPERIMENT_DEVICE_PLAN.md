# Acoustic single-Rayleigh BIC: PRL figure, experiment, and device plan

## 1. Scope and fixed design

All final scattering and field predictions in the figure gallery use the
manufacturable 200 kHz design

\[
a=6.871135\ \mathrm{mm},\quad
\theta_{\mathrm{BIC}}=5.251216^\circ,
\]

\[
(d_1,d_2)=(4.541418,1.080580)\ \mathrm{mm},
\]

\[
(w_1,w_2)=(3.964770,1.053048)\ \mathrm{mm},\quad
g=0.926658\ \mathrm{mm}.
\]

The strict homogeneous mode is evaluated at \(N=313,K=39\). The independent
pole continuation uses the exact moderate root \(N=121,K=15\), and is shown as
a convergence diagnostic rather than being silently mixed with driven
scattering data.

## 2. How the acoustic figures follow arXiv:2607.07228

The arXiv paper uses four steps:

1. meta-atom response, angle-frequency map, a fixed-angle cut, and the number
   of open channels;
2. real pole dispersion, imaginary linewidth/Q, and fixed-angle spectra;
3. Rayleigh lines, open-channel topology, and robustness to added radiation
   degrees of freedom;
4. a full-wave map, Q curve, and representative eigenfields.

The acoustic translation is:

| arXiv logic | Acoustic observable | Gallery figure |
|---|---|---|
| single-particle/anapole response | two-groove radiation cancellation and small-cavity phase control | Fig. 7 |
| reflectance map | \(A_n\), \(\eta_n\), reflection phase, and groove-field maps | Fig. 4 |
| fixed-angle reflectance cuts | \(\arg A_0\) and \(\eta_{-1}\) cuts | Fig. 3(c,d) |
| open-channel count | \(n=-1,0,+1\) propagating/Rayleigh/evanescent topology | Fig. 1(b,c) |
| pole reaches Rayleigh anomaly | \(\operatorname{Re}\Omega_p\) versus \(\kappa\) | Fig. 3(a) |
| linewidth collapse and Q divergence | \(\operatorname{Im}\Omega_p\) and \(Q\) | Fig. 3(b) |
| full-wave field | complex pressure and time-averaged intensity | Fig. 5 |

The optical paper's cubic Q law must not be copied into the acoustic claim.
The present finite modal model shows an approximately quadratic trend; a
universal exponent has not been established.

## 3. Ten multi-panel article figures

1. **Structure and single-Rayleigh topology.** Physical unit cell, Rayleigh
   lines, channel status, and the \(\pm\theta\) time-reversal pair.
2. **Strict BIC evidence.** Root convergence, strict singular residual,
   \(A_{-1}=A_0=0\), and exterior square-integrability.
3. **Pole and scattering cuts.** Pole-to-branch-point trajectory, Q collapse,
   specular phase cuts, and anomalous-order cuts.
4. **Theoretical water-tank observable atlas.** Maps of
   \(|A_{-1}|,|A_0|,|A_{+1}|,\eta_{-1},\eta_0\), and groove-field norm.
5. **Fields.** Off-resonance, ideal high-Q quasi-BIC, practical 30% router,
   and dark homogeneous BIC eigenfield.
6. **Reciprocal signed-order routing.** Positive/negative angle pair, routing
   efficiency, full power-normalized scattering matrices, reciprocity and
   unitarity residuals.
7. **Why the small cavity is essential.** Single-/two-cavity comparison,
   truncation residuals, and separate complex radiation phasors for
   \(A_0\) and \(A_{-1}\).
8. **Manufacturability and robustness.** Physical dimensions, tolerance map,
   sound-speed drift, and retuned device-point convergence.
9. **Water-tank experiment.** Tank geometry, coarse-to-fine scan, complex
   field-to-Floquet extraction, and the full evidence chain.
10. **Device demonstration.** A theta-frequency regime map separates the
    strict BIC, the near-Rayleigh 30.6% signed-order router, and a distinct
    same-cell off-Rayleigh point with near-unity anomalous reflection.
    Modal-truncation convergence, frequency/angle windows, and reciprocal
    output directions make the device claim independently testable.

For a five-figure PRL main text, the recommended sequence is Fig. 1, Fig. 2,
Fig. 3, Fig. 5, and a combined experimental version of Fig. 6/9/10. Figures
4, 7, and 8 naturally belong in Supplemental Material unless the small-cavity
mechanism becomes a main novelty.

## 4. Water-tank experiment

### 4.1 Samples and controls

- Full two-groove arrays with \(M=20,40,60\) periods.
- A flat rigid reference plate.
- A single-wide-groove sample.
- A full sample with the small groove blocked.
- Optional small-groove width/depth perturbation samples.

The widths of the \(M=40\) and \(M=60\) arrays are about 274.85 mm and
412.27 mm. The sample material and backing thickness must be checked with a
fluid-structure model; “rigid” should not be assumed solely from the material
name.

### 4.2 Calibration

1. Measure water temperature and the actual sound speed by time of flight.
2. Measure the transmitter's complex free-field pattern and angular spectrum.
3. Measure the hydrophone transfer function.
4. Measure the flat-plate complex reflection reference.
5. Time-gate tank-wall and surface echoes.

If the measured sound speed is 1480 m/s rather than 1500 m/s, the 200 kHz
design feature moves to approximately 197.33 kHz.

### 4.3 Scan sequence

1. Coarse map: \(190\text{--}205\) kHz and
   \(-7^\circ\leq\theta\leq7^\circ\), with 0.1--0.2 kHz frequency steps.
2. Fine map: around the measured Rayleigh/quasi-BIC feature, use 10--50 Hz
   steps and a fine angle scan around \(\pm5.25^\circ\) and \(\pm6.0^\circ\).
3. Independent off-Rayleigh device map: scan
   \(202.30\text{--}202.56\) kHz for
   \(25^\circ\leq|\theta|\leq40^\circ\), then refine around
   \(|\theta|=32.6328^\circ\) with 10--50 Hz steps. The expected target order
   is \(n=-1\) for \(+\theta\) and \(n=+1\) for \(-\theta\), with output
   angles near \(-32.63^\circ\) and \(+32.63^\circ\), respectively.
4. Near field: scan complex pressure over the central 3--5 periods, from
   approximately 0.25 mm to 5--10 mm above the surface, with 0.25--0.5 mm
   lateral sampling where possible.
5. Far field: scan an angular arc containing the specular beam and the
   near-grazing anomalous beam; for the off-Rayleigh point, include the
   reciprocal output lobes near \(\theta_o=\mp32.63^\circ\).
6. Ring-down: excite with a long narrowband burst, switch off, and fit the
   field amplitude \(e^{-t/\tau_A}\), using \(Q=\pi f\tau_A\).

### 4.4 Extracting diffraction amplitudes

Measure the incident field without the sample and the total field with the
sample, then form the complex scattered field. At a fixed height \(y_0\), fit

\[
p_{\mathrm{sc}}(x,y_0)=\sum_n\widetilde A_n(y_0)e^{ik_{x,n}x}
\]

over the central aperture. Propagate the fitted coefficients back to the
surface and calculate

\[
\eta_n=\frac{\operatorname{Re}k_{y,n}}{k_{y,\mathrm{inc}}}|A_n|^2.
\]

The total driven-scattering \(A_0\) contains a nonresonant rigid-background
reflection. It is not the homogeneous BIC coefficient \(A_0^{\mathrm{BIC}}=0\).
For experimental pole evidence, fit

\[
A_n(f)=A_{n,\mathrm{bg}}(f)+\frac{r_n}{f-f_p}
\]

and track the pole/residue rather than equating a total reflection minimum
with the BIC radiation zero.

## 5. Solid device conclusions

### 5.1 Current robust result

At the practical high-truncation operating point

\[
\theta=\pm6.000^\circ,\qquad f=200.216874\ \mathrm{kHz},
\]

the positive-angle input gives

\[
\eta_{-1}=0.30613,\qquad \eta_0=0.69387,
\]

with the anomalous beam at approximately \(-80.34^\circ\). Reversing the
incidence angle gives the reciprocal partner

\[
\eta_{+1}=0.30613,
\]

at approximately \(+80.34^\circ\).

After separately retuning the operating frequency at
\(K=15,23,35,39\), the anomalous efficiency converges as

\[
0.30737,\ 0.30662,\ 0.30620,\ 0.30613,
\]

while the frequency converges from 200.21475 to 200.21687 kHz. This remains the
strongest current **Rayleigh-linked** device number that is sufficiently
converged to use as a theoretical prediction; the separate off-Rayleigh
near-unity point is documented next.

### 5.2 Same-cell near-unity off-Rayleigh reflection

The same two-groove geometry also has a separate, high-efficiency operating
point at

\[
\theta_i=\pm32.6328^\circ,\qquad f=202.430\ \mathrm{kHz}.
\]

For \(+\theta_i\), the selected order is \(n=-1\) and the output is
\(\theta_o\approx-32.63^\circ\); for \(-\theta_i\), reciprocity maps it to
\(n=+1\) and \(\theta_o\approx+32.63^\circ\). At \(N=401,K=50\), the modal
calculation gives

\[
\eta_{\mathrm{target}}=0.9996185,\qquad
\eta_0=0.0003815,\qquad
\delta_{\mathrm{energy}}=4.4\times10^{-16},
\]

with condition number \(7.3\times10^3\). The convergence sequence for
\((N,K)=(81,11),(121,15),(161,20),(201,25),(313,39),(401,50)\) is
\(0.97791,0.99298,0.99864,0.99963,0.99992,0.99962\). In the plotted
\(N=201,K=25\) frequency cut, the \(\eta>0.90\) window is
\(202.365\text{--}202.491\) kHz and the \(\eta>0.99\) window is
\(202.408\text{--}202.445\) kHz. At 202.430 kHz, the corresponding
\(N=201,K=25\) angle windows are approximately
\(27.6^\circ\text{--}38.0^\circ\) and
\(30.7^\circ\text{--}34.6^\circ\), respectively.

This point is **off Rayleigh**: the active \(n=-1\) threshold at
\(|\theta|=32.6328^\circ\) is approximately 141.825 kHz, so the operating
frequency is about 60.6 kHz above it. Therefore the 99.96% result must not be
called BIC-enhanced, Rayleigh-critical, or nonreciprocal. It is an independent
two-channel anomalous-reflection state of the same manufacturable unit cell.

### 5.3 Ideal limit versus experiment

A point closer to the BIC, near

\[
\theta=5.27044^\circ,\qquad f=200.00578\ \mathrm{kHz},
\]

gives about 85.7% anomalous conversion in the lossless high-truncation model.
However, its matrix condition number is approximately \(1.6\times10^{10}\),
its field enhancement is enormous, and its usable bandwidth tends to zero as
the exact BIC is approached. It is an ideal-limit result, not yet an
experiment-ready claim. Viscous loss, elastic compliance, finite aperture,
and fabrication disorder must be included before promoting this number.

### 5.4 Reciprocity no-go

For the two open channels, use the power-normalized scattering matrices. Time
reversal gives

\[
S_+(\kappa)=S_-^T(-\kappa),
\]

and the present lossless backed structure also obeys

\[
S_\pm^\dagger S_\pm=I.
\]

Consequently,

\[
\eta_{-1}(+\theta)=\eta_{+1}(-\theta).
\]

The system is therefore a **reciprocal signed-order router**, not a
nonreciprocal isolator. “One anomalous order for a fixed input sign” is valid;
“one-way propagation” is not.

### 5.5 Routes beyond the present device

- **Higher reciprocal conversion:** the off-Rayleigh point above already
  supplies a near-unity two-channel benchmark. A separate route is to jointly
  optimize a deliberately broadened quasi-BIC and specular cancellation. Unit
  conversion is allowed, but it occurs for both time-reversed directions.
- **Directional but reciprocal response:** introduce controlled loss or a
  third radiative port. This can make the two same-order incidence experiments
  look asymmetric, while the exact time-reversed scattering process remains
  reciprocal.
- **True isolation:** add time modulation, mean flow, rotation, or another
  genuine time-reversal-breaking bias. This is a new physical platform rather
  than a claim supported by the present passive grating.

## 6. Files

- Figure generator: `Ni2019_MATLAB/run_PRL_candidate_figure_gallery.m`
- Route objective: `Ni2019_MATLAB/ni2019_route_objective.m`
- Cached data: `Ni2019_MATLAB/results/PRL_article_figure_data_v5.mat`
- Gallery: `Ni2019_MATLAB/results/PRL_gallery/`
- Contact sheet: `Ni2019_MATLAB/results/PRL_gallery/PRL_Figure_gallery_contact_sheet.png`
