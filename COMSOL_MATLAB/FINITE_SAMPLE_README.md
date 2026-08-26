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
- A centered `80 mm` active line imposes a phased Gaussian normal acceleration
  with a `25 mm` waist. The apodization suppresses the sinc sidelobes of a
  hard-edged aperture; it does not remove the finite beam's angular spread.
  The remaining exterior boundaries use plane-wave radiation conditions.

## Reproduction

Run, in this order, with MATLAB LiveLink for COMSOL:

1. `build_finite_20period_80mm_source.m`
2. `scan_finite_20period_theta_frequency.m`
3. `export_finite_20period_bic_fields.m`
4. `build_finite_incident_reference.m`
5. `export_finite_gaussian_three_angle_fields.m`
6. `../Ni2019_MATLAB/run_fig4_wideangle_pole_180k.m`
7. `plot_comsol_numerical_experiment.py`

The first script creates the complete finite model and a coarse frequency
sweep. The second resolves five angular cuts with a 2-Hz frequency step. The
third exports fields at `f0-5 Hz`, `f0`, and `f0+5 Hz`, where
`f0 = 180.046681327 kHz`. The Python script creates the vector PDF/SVG and the
600-dpi preview used for the numerical-experiment figure.

The incident-reference model removes the sample, keeps the same Gaussian
source, and replaces the lower surface by a radiation boundary. The plotted
reflected field is calculated without visual post-processing as
`p_reflected = p_total - p_incident` on the common COMSOL grid. The three
representative states follow the physical eigenpole branch: `0 deg` at
`178.7802 kHz`, the rounded-cell Rayleigh point `7.05 deg` at
`180.046681 kHz`, and `10.054 deg` at `180.498322 kHz`.

The wide-angle continuation covers approximately `0--20 deg`. Its computed
radiative linewidth is zero at the strict Rayleigh BIC (`7.0996629 deg`), about
`1.15 Hz` at `8.00 deg`, `15.46 Hz` at `10.05 deg`, and `147.46 Hz` at
`20.03 deg`. The same eigenbranch also reaches a distinct single-open-channel
dark state at the Gamma point. Because the unit cell is not mirror symmetric,
this endpoint is not labelled as a symmetry-protected BIC.

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

The field panels therefore are not interpreted by their absolute darkness.
They expose the incident/reflected separation and provide the complex scan-line
field used for an angular-spectrum decomposition. The direct scattering
signature is the linewidth of the background-subtracted Floquet projection:
its FWHM agrees with the outgoing eigenpole linewidth and approaches zero at
the Rayleigh point. A nonzero background reflection at the BIC is expected and
must not be misidentified as radiative leakage of the bound state.

The large `.mph` files are published as release assets at:

<https://github.com/YongLi2010/UnderwaterRayleighBIC/releases/tag/comsol-finite-20p-v1>
