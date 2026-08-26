# Figure QA notes

Core conclusion: the outgoing eigenpole and the linewidth extracted from
driven reflected fields collapse at the same Rayleigh point, while a finite
20-period Gaussian-beam COMSOL model visualizes how the reflected far field
changes as the second Floquet channel opens.

| Panel | Unique claim | Source | Visual audit |
|---|---|---|---|
| a | Finite sample, Gaussian aperture, and scan geometry | Model parameters | Pass; no collisions |
| b | The resonance linewidth visibly pinches off at the Rayleigh BIC | `wideangle_physical_pole.csv` | Pass; pole-derived Lorentzian spectral function is explicitly identified |
| c | Quantitative linewidth evolution from 0 to 20 degrees | `wideangle_physical_pole.csv` | Pass; logarithmic axis and zero-value display floor documented in source |
| d | Incident and reflected fields at normal incidence, the Rayleigh point, and above it | three finite-sample COMSOL solutions and matched no-sample references | Pass; common grid and one physical color scale |
| e | Channel opening and scattering-derived linewidth | scan-line angular spectra and background-subtracted driven Floquet projections | Pass; propagating interval and independent eigenpole curve are shown |

Exports: editable SVG, editable PDF, 600-dpi PNG preview, and 600-dpi TIFF.
The PDF text audit passes with a minimum detected glyph size of `5.25 pt`.

The static validator's remaining warnings are reviewed false positives:

- mathematical legend labels do not require sentence-style capitalization;
- the plotted logarithm is explicitly guarded by `maximum(value, 1e-30)`;
- the rendered PDF glyph audit, rather than the conservative source estimate,
  verifies that all glyphs exceed 5 pt.

No experimental data or stochastic uncertainty is presented. Panel b is a
pole-derived normalized spectral function; panels a and d use deterministic
COMSOL geometry and fields. Panel c uses the directly computed complex
eigenfrequency linewidth. In panel e, the orange symbols are theory-derived
FWHM values from complex driven-scattering amplitudes, not experimental data.

Image integrity: the field maps show the real pressure without local contrast
adjustment, stitching, or selective masking. All six maps share one color
scale. The reflected field is obtained by complex subtraction on a common
grid. The angular spectra use one Hann window over the full 148.4-mm scan line
and are normalized only to each spectrum's own maximum for channel visibility.

## Figure 5: ideal-versus-finite scattering

The companion comparison figure separates three effects that were previously
visually conflated:

| Panel | Unique claim | Source | Visual audit |
|---|---|---|---|
| a | The ideal infinite system produces pure Bloch-reflected fields | full modal-matching reflected fields | Pass; five displayed periods, no local normalization |
| b | Finite aperture and edge scattering broaden the reflected field | three COMSOL sample/reference pairs | Pass; one incident normalization and one color scale |
| c | Ideal Floquet lines observed through a finite measurement window | complex modal-matching Floquet coefficients, reconstructed at 20 mm | Pass; the same 20-period Hann-window projection as panel d is used |
| d | A finite sample and Gaussian beam broaden the angular spectrum | dense windowed Fourier integral of COMSOL scan-line fields | Pass; the same axes, window length, and normalization as panel c are used |
| e | The ideal plane-wave resonance pinches off at the BIC | pole-derived spectral density | Pass; logarithm is guarded |
| f | The 25-mm-waist beam averages nearby leaky Bloch angles | Gaussian angular convolution of the same spectral density | Pass; same axes and color scale as panel e |
| g | Finite angular width regularizes the vanishing linewidth | eigenpole, ideal driven-scattering FWHM, and Gaussian-convolved FWHM | Pass; logarithmic axis spans all three scales |

The first two infinite-period columns use a `+0.1 Hz` limiting-background
solution because the exact forced problem is non-unique at a dark state. This
regularization is stated in the reproduction notes and is not presented as a
finite BIC linewidth. The displayed total background field is distinct from
the homogeneous BIC eigenfield.

Panel d is extracted directly from the complex finite-sample COMSOL scan-line
fields. Panels f and g are a predicted finite-beam frequency scan: the ideal
pole-derived spectral density is averaged with the angular intensity of the
`exp[-(x/w)^2]`, `w = 25 mm` source. Thus the approximately `600 Hz` apparent
width is an aperture-averaging prediction, not a fitted linewidth from the
three static COMSOL snapshots and not an experimental observation. The
80-mm source truncation is neglected in this convolution because its Gaussian
amplitude is only about `0.077` at each edge.

Panels c and d now use an identical observation operator. For panel c, the
ideal infinite-periodic pressure is reconstructed from the complex Floquet
coefficients at `y = 20 mm`, sampled across `20a`, multiplied by a Hann window,
and projected on a dense `kx/k0` grid. Panel d applies that projection directly
to the finite-sample COMSOL reflected field. Thus the width of the ideal peaks
in panel c is solely the selected observation-window resolution; in the
infinite-window limit those peaks become discrete Floquet delta lines. The
plotted data for both panels are exported in `Fig5_cd_angular_spectra.csv`.
