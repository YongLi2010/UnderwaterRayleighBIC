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

## Finite-versus-infinite comparison

The companion comparison figure separates three effects that were previously
visually conflated:

| Panel | Unique claim | Source | Visual audit |
|---|---|---|---|
| a | Finite aperture and edge scattering broaden the measured field | three COMSOL sample/reference pairs | Pass; one incident normalization and one color scale |
| b | An infinite Bloch plane wave has no aperture diffraction | full modal-matching reflected fields | Pass; five displayed periods, no local normalization |
| c | A finite scan produces a continuous but resolution-limited angular spectrum | dense windowed Fourier integral of COMSOL scan-line fields | Pass; no raw FFT-bin joining |
| d | The infinite far field consists of discrete Floquet lines | exact modal-matching power fractions | Pass; stems are used instead of continuous interpolation |
| e | Channel opening alone is not the BIC criterion | pole linewidth and driven-scattering FWHM | Pass; independent quantities agree |

The first two infinite-period columns use a `+0.1 Hz` limiting-background
solution because the exact forced problem is non-unique at a dark state. This
regularization is stated in the reproduction notes and is not presented as a
finite BIC linewidth. The displayed total background field is distinct from
the homogeneous BIC eigenfield.
