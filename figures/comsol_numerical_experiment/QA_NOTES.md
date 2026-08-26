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
