# Figure QA notes

Core conclusion: the outgoing eigenpole and the linewidth extracted from
driven reflected fields collapse at the same Rayleigh point, while a finite
20-period Gaussian-beam COMSOL model visualizes how the reflected far field
changes as the second Floquet channel opens.

| Panel | Unique claim | Source | Visual audit |
|---|---|---|---|
| a | Ideal Bloch-wave drive and Fourier projection onto the reflected orders | periodic modal model | Pass; no finite-aperture assumptions |
| b | The resonance linewidth visibly pinches off at the Rayleigh BIC | `wideangle_physical_pole.csv` | Pass; pole-derived Lorentzian spectral function is explicitly identified |
| c | Quantitative linewidth evolution from 0 to 20 degrees | `wideangle_physical_pole.csv` | Pass; logarithmic axis and zero-value display floor documented in source |
| d | Extracted pole linewidth and integrated weights of the two reflected orders | adaptive periodic scattering response | Pass; both channels are retained |
| e | Near-field scan and Floquet reconstruction at a nearby driven state | periodic reference-subtracted field | Pass; common color scale and explicit channel markers |

Exports: editable SVG, editable PDF, 600-dpi PNG preview, and 600-dpi TIFF.
The PDF text audit passes with a minimum detected glyph size of `5.25 pt`.

The static validator's remaining warnings are reviewed false positives:

- mathematical legend labels do not require sentence-style capitalization;
- the plotted logarithm is explicitly guarded by `maximum(value, 1e-30)`;
- the rendered PDF glyph audit, rather than the conservative source estimate,
  verifies that all glyphs exceed 5 pt.

No experimental data or stochastic uncertainty is presented. Panel b is a
pole-derived normalized spectral function; panel a is a schematic readout of
the ideal periodic problem. Panels c and d use directly computed complex
driven-scattering amplitudes, and panel e uses a deterministic reference-
subtracted field. All are theory predictions, not experimental data.

Image integrity: the field maps show the real pressure without local contrast
adjustment, stitching, or selective masking. All six maps share one color
scale. The reflected field is obtained by complex subtraction on a common
grid. The angular spectra use one Hann window over the full 148.4-mm scan line
and are normalized only to each spectrum's own maximum for channel visibility.

## Figure 5: Gaussian-beam regularization of the Rayleigh BIC

Core conclusion: a finite Gaussian beam does not directly reproduce the
single-Bloch-wave BIC limit, but its calibrated angular spectrum supplies a
controlled route to recover the two channel responses and extrapolate toward
that limit.

| Panel | Unique claim | Source | Visual audit |
|---|---|---|---|
| a | Complex reflected-field scanning can recover `r0(kx,f)` and `r-1(kx,f)` | current 20-period, 80-mm-aperture measurement geometry | Pass; schematic dimensions match the model |
| b | The present `w=25 mm` source has a broad `7.2 deg` intensity FWHM | analytic Fourier transform of `exp[-(x/w)^2]` | Pass; all curves use the same normalization |
| c | A single Bloch wave resolves the pole-linewidth pinch at the BIC | outgoing-eigenpole frequency and linewidth | Pass; guarded logarithm and exact BIC marker |
| d | The current Gaussian beam averages the pinch over neighboring leaky states | angular convolution of the same pole spectral density | Pass; axes and color scale match panel c |
| e | Angular averaging creates a finite spectral-resolution floor | pole linewidth and Gaussian-weighted frequency-distribution width | Pass; the floor is not labelled as an intrinsic linewidth |
| f | Increasing waist systematically lowers that floor | the same forward model over `w=25--300 mm` | Pass; current aperture is explicitly identified |

The active source is `80 mm` wide, whereas the Gaussian waist is `25 mm`.
For the pressure amplitude `exp[-(x/w)^2]`, the `1/e^2` intensity diameter is
`2w=50 mm`; at the source edges the amplitude is approximately `0.077`.

Panels c and d use a shared absolute spectral-density scale. Panel d is a
forward prediction, not COMSOL data or experimental data. The finite-beam
spectral width in panels e and f is the Gaussian-equivalent FWHM of the pole
center-frequency distribution induced by the known angular spectrum. This
moment-based quantity is used because the fully convolved response can become
multipeaked, making an envelope FWHM discontinuous and estimator-dependent.
It therefore represents an angular-dispersion resolution floor, not a fitted
resonance linewidth.

The wider-waist curves assume proportionally enlarged, weakly truncated
sources. Keeping the present edge amplitude requires an active width of about
`3.2w`; the sample must also be enlarged accordingly. They are design
predictions and must not be described as results from the current 80-mm source.
Quantitative source data are exported in `Fig5_gaussian_beam_scaling.csv` and
`Fig5_minimum_width_vs_waist.csv`.

Exports: editable SVG and PDF plus 600-dpi PNG/TIFF. The final PDF text audit
passes with a minimum detected glyph size of `5.04 pt`. The static logarithm
warning is a reviewed false positive: both spectral maps explicitly apply a
`1e-6` positive display floor before `log10`.
