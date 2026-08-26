# Finite 20-period COMSOL numerical experiment

This model mimics a two-dimensional underwater reflection experiment for the
rounded, fabrication-ready Rayleigh-BIC sample.

## Geometry and excitation

- 20 periods, `a = 7.42 mm`; total sample width `148.40 mm`.
- Large groove: `w1 = 4.41 mm`, `d1 = 5.15 mm`.
- Small groove: `w2 = 1.00 mm`, `d2 = 1.32 mm`.
- Intra-cell and inter-cell gaps: `g12 = 1.00 mm`, `g21 = 1.01 mm`.
- Stainless-steel plate thickness: `10.00 mm`.
- Water: `c = 1500 m/s`, `rho = 1000 kg/m^3`.
- Two-dimensional model; the grooves and plate are invariant along `z`.
- The steel is represented geometrically, while the acoustic calculation uses
  the sound-hard limit at the water--steel interface.
- A centered `80 mm` active line imposes a phased normal acceleration. The
  remaining exterior boundaries use plane-wave radiation conditions.

## Reproduction

Run, in this order, with MATLAB LiveLink for COMSOL:

1. `build_finite_20period_80mm_source.m`
2. `scan_finite_20period_theta_frequency.m`
3. `export_finite_20period_bic_fields.m`
4. `plot_comsol_numerical_experiment.py`

The first script creates the complete finite model and a coarse frequency
sweep. The second resolves five angular cuts with a 2-Hz frequency step. The
third exports fields at `f0-5 Hz`, `f0`, and `f0+5 Hz`, where
`f0 = 180.046681327 kHz`. The Python script creates the vector PDF/SVG and the
600-dpi preview used for the numerical-experiment figure.

## Physical interpretation

The infinite-period port model resolves an extremely narrow quasi-BIC feature:
at `theta = 7.10 deg`, the sampled cavity-pressure peak has a width no larger
than the `0.5 Hz` frequency step. The rounded Rayleigh crossing associated with
`f0` occurs near `theta = 7.054 deg`.

The finite 20-period model with an 80-mm aperture does **not** reproduce an
infinite-Q singularity. Its aperture is only about ten wavelengths wide, so the
incident beam contains a broad angular spectrum, and both the finite length and
finite aperture regularize the response. This is the experimentally meaningful
limit: report a minimum resolvable linewidth, not a directly measured
mathematical BIC.

The large `.mph` files are published as release assets at:

<https://github.com/YongLi2010/UnderwaterRayleighBIC/releases/tag/comsol-finite-20p-v1>
