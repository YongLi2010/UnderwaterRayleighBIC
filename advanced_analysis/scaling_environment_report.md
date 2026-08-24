# Scaling and sound-speed robustness of the strict Rayleigh BIC

This analysis supplements the existing theoretical manuscript without changing its physical claim. It uses the existing pole-free homogeneous modal-matching operators in `Ni2019_MATLAB/ni2019_strict_rayleigh_operator.m` and `Ni2019_MATLAB/ni2019_full_eigen_operator.m`. No fitted radiation model or surrogate solver is used.

The reproducible entry point is:

```text
matlab -batch "run('advanced_analysis/scaling_environment.m')"
```

The script reads the full-precision retained 200-kHz, 1-mm-constrained design from `Ni2019_MATLAB/results/StrictRayleighBIC_200kHz_min1mm.mat`, evaluates the nominal point, performs the two scaling sweeps, and writes all output to `advanced_analysis/scaling_environment/`.

## 1. Normalization and exact similarity law

The solver uses

\[
 \Omega=\frac{af}{c},\qquad
 \kappa=\Omega\sin\theta,
 \qquad
 \frac{k_{y,n}}{k_0}=\sqrt{1-\left(\sin\theta+\frac{n}{\Omega}\right)^2},
\]

with the outgoing square-root branch. For positive incidence at the nominal point, the active `n=-1` Rayleigh condition is

\[
 \Omega=1-\kappa,
 \qquad \kappa+\Omega=1.
\]

The nominal normalized values are:

```text
Omega0   = 0.916151303447099
kappa0   = 0.0838486965529006
theta0   = 5.251216456 deg
d/a      = [0.660941508245, 0.157263731977]
w/a      = [0.577018249316, 0.153256794269]
g/a      = 0.134862478208
```

The dimensional nominal values are:

```text
c0       = 1500 m/s
f0       = 200.000000 kHz
a0       = 6.871134776 mm
lambda0  = 7.500000000 mm
```

If all normalized parameters are held fixed, let

\[
 a=s a_0,\qquad
 (d,w,g)=s(d_0,w_0,g_0),\qquad
 \theta=\theta_0.
\]

Then the Rayleigh-BIC frequency is exactly

\[
 f_{\mathrm{BIC}}(c,s)=\frac{\Omega_0c}{s a_0}
 = f_0\frac{c}{c_0}\frac{1}{s}.
\]

The code nevertheless rebuilds and evaluates the full modal operator at each dimensional point. Thus the linear and inverse-linear trends shown in the figures are checked against the actual modal matching calculation rather than inserted as plotted formulas.

## 2. Nominal reproduction

At `(N,K)=(313,39)`, the script reproduces the strict homogeneous operator at the final retained geometry:

| quantity | value |
|---|---:|
| target threshold order | `n=-1` |
| removed open/threshold orders | `[-1, 0]` |
| `sigma_min/sigma_max` | `1.1817e-15` |
| reduced raw residual | `9.6944e-14` |
| nominal threshold status | true |
| finite-flux open power | zero to numerical precision |

The small difference from the normalized-coordinate value quoted in the manuscript (`1.2604e-15` and `1.3766e-13`) comes from rebuilding the same matrix in physical length units. Both calculations use the full-precision geometry and the same truncation; the machine-level null and exact removed-channel cancellation are unchanged.

## 3. Sound-speed scaling

For `s=1`, the sweep `c=1400:10:1600 m/s` retunes the frequency according to `f_BIC=Omega0*c/a0`, while keeping the physical geometry and incidence angle at their nominal values. The resulting range is

```text
f_BIC(1400 m/s) = 186.666667 kHz
f_BIC(1500 m/s) = 200.000000 kHz
f_BIC(1600 m/s) = 213.333333 kHz
```

Across all 21 points:

```text
strict raw residual range = 9.6944e-14 ... 1.1388e-13
singular-value ratio      = 2.1252e-16 ... 1.1817e-15
max |Omega-Omega0|        = 1.11e-16
```

The residual variation is numerical matrix-conditioning/truncation noise, not a systematic loss of the cancellation. The relevant conclusion is conditional: within the lossless, nondispersive, infinite-period model, a similarity-retuned design preserves the strict Rayleigh-BIC mechanism under a change of sound speed.

The generated data are in `scaling_environment/fBIC_vs_sound_speed.csv`.

## 4. Geometric scale sweep

At `c=c0`, the scale factor was varied over `s=0.50:0.10:2.00`, with every physical length scaled by `s` and the frequency retuned by `1/s`:

```text
f_BIC(s=0.50) = 400.000000 kHz
f_BIC(s=1.00) = 200.000000 kHz
f_BIC(s=2.00) = 100.000000 kHz
```

Over the 16 evaluated scale factors:

```text
strict raw residual range = 9.1617e-14 ... 1.8606e-13
singular-value ratio      = 1.7131e-16 ... 1.8143e-15
```

The result is the expected dimensional similarity: the normalized matrix, channel topology, and cancellation constraints are unchanged when all lengths and the wavelength are scaled together. It is not a claim that losses, finite apertures, transducer coupling, or fabrication errors are scale independent.

The generated data are in `scaling_environment/fBIC_vs_scale_factor.csv`.

## 5. Fixed-device sound-speed control

As a control, the actual geometry, `f=200 kHz`, and `theta=theta0` were held fixed while `c` was swept. Now

\[
 \Omega_{\mathrm{fixed}}(c)=\frac{a_0f_0}{c},
 \qquad
 \kappa_{\mathrm{fixed}}(c)=\Omega_{\mathrm{fixed}}\sin\theta_0,
\]

so the nominal `n=-1` threshold is met only at `c=1500 m/s`. Away from that point the target order is either a finite-flux propagating channel (`c<1500` m/s) or an evanescent order (`c>1500` m/s); it is not legitimate to call the result a Rayleigh BIC.

The full, unreduced homogeneous operator was therefore used for this control. Its smallest-singular-vector mode gives the ordinary real-frequency operator residual and the normalized open-channel amplitude. Selected values are:

| `c` (m/s) | `Omega` | Rayleigh detuning `||kappa-1|-Omega|` | full raw residual | open-amplitude ratio |
|---:|---:|---:|---:|---:|
| 1400 | 0.98159068 | 0.07142857 | 0.60409 | 0.77557 |
| 1500 | 0.91615130 | `0` | `8.80e-14` | `4.54e-16` |
| 1600 | 0.85889185 | 0.06250000 | 0.36823 | 0.09491 |

Across the complete fixed-device sweep, the full-operator raw residual ranges from `8.80e-14` to `6.49e-1`, and the open-amplitude ratio ranges from `4.54e-16` to `7.76e-1`. These values show why the exact retuning statement must be separated from a fixed-frequency environmental tolerance statement.

The script also records a `forced_zero_control_residual`. This is obtained by applying the existing strict operator while deleting the nominal target amplitude even when the target is no longer grazing. It is retained only as an over-constrained numerical control and must not be interpreted as a BIC criterion away from a Rayleigh point.

The generated data are in `scaling_environment/fixed_device_sound_speed_control.csv`.

## 6. Figures and data products

The MATLAB script writes vector PDF figures using `exportgraphics(..., 'ContentType','vector')`:

```text
scaling_environment/fig_scaling.pdf
scaling_environment/fig_environment.pdf
```

`fig_scaling.pdf` contains `f_BIC` versus sound speed and geometric scale factor. `fig_environment.pdf` contains the retuned frequency, retuned homogeneous residual, fixed-device residual/open amplitude, and the crossing of `Omega` with the Rayleigh threshold.

Machine-readable outputs are:

```text
scaling_environment/fBIC_vs_sound_speed.csv
scaling_environment/fBIC_vs_scale_factor.csv
scaling_environment/fixed_device_sound_speed_control.csv
scaling_environment/scaling_environment_data.mat
```

The PDF files were rendered and visually checked after generation; Poppler reports one page per figure and no raster image objects, consistent with vector output.

## 7. Manuscript interpretation and limitations

The safe manuscript-level statement is that the strict Rayleigh-BIC calculation is dimensionally scale invariant under exact similarity retuning, and that a fixed device has a frequency/angle-dependent tolerance to sound-speed changes because the normalized Rayleigh condition moves.

This calculation does not establish scale invariance of a real water-tank experiment. It omits thermoviscous loss, material dispersion, fluid-structure compliance, finite sample width, source/aperture effects, and calibration errors. Those effects should be included in a later finite-device study.
