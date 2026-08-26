# Figure QA notes

Core conclusion: the periodic COMSOL model exhibits a collapsing scattering
linewidth near the Rayleigh crossing, whereas the finite 20-period, 80-mm
aperture model provides the regularized field and hydrophone observables of a
numerical experiment.

| Panel | Unique claim | Source | Visual audit |
|---|---|---|---|
| a | Finite sample, aperture, and scan geometry | Model parameters | Pass; no collisions |
| b | Global resonance trajectory and Rayleigh line | `comsol_coarse_angle_frequency_map.csv` | Pass; monotone axes checked |
| c | Sub-hertz channel conversion and cavity enhancement | `comsol_fine_frequency_cuts.csv` | Pass; no interpolation |
| d | Finite-sample total pressure field | state-2 real/imaginary grids | Pass; symmetric clipping at 99.2 percentile documented in source |
| e | Directly measurable hydrophone scan | `finite_bic_state_2_scanline.csv` | Pass; source footprint shown explicitly |

Exports: editable SVG, editable PDF, 600-dpi PNG preview, and 600-dpi TIFF.
The PDF text audit passes with a minimum detected glyph size of `5.25 pt`.

The static validator's remaining warnings are reviewed false positives:

- mathematical legend labels do not require sentence-style capitalization;
- the plotted logarithm is explicitly guarded by `maximum(value, 1e-30)`;
- the rendered PDF glyph audit, rather than the conservative source estimate,
  verifies that all glyphs exceed 5 pt.

No experimental data or stochastic uncertainty is presented. All plotted
quantities are deterministic COMSOL outputs.
