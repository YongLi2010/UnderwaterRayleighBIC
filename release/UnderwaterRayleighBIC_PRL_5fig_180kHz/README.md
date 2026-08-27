# Underwater Rayleigh BIC — five-figure PRL resource package

This package contains the current, internally consistent 180-kHz design chain
(`a = 7.4167 mm`, `c = 1500 m s^-1`, `theta_BIC = 7.09966 deg`).  The paper
source is `main.tex`; the detailed derivations and numerical checks are in
`suppl.tex`; references are in `ref.bib`.

The vector figures are supplied as both PDF and editable SVG:

1. `fig1_potential_refined_v2`: conventional BIC, Rayleigh anomaly, and
   Rayleigh BIC concept.
2. `fig2_rayleigh_bic_eigenproblem`: pole spectrum, Rayleigh crossing,
   linewidth divergence, channel zeros, and homogeneous fields.
3. `fig3_common_dark_state_nature`: two-mode/two-channel CMT mechanism.
4. `fig4_scattering_observability_nature`: ideal Bloch-wave scattering and
   channel-resolved observables.
5. `fig5_finite_beam_observability`: finite Gaussian-beam angular averaging
   and the resulting resolution floor.

`main.tex` currently uses Figures 1–4 as the compact Letter core.  Figure 5
is included as an optional theory figure for the finite-beam measurement
discussion; it is deliberately not presented as experimental data.

The separate 200-kHz root is not mixed into this package.  All plotted values
are theory predictions; no experimental measurements are claimed.

To compile from this directory, use the standard RevTeX/BibTeX sequence:

```sh
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
pdflatex suppl.tex && bibtex suppl && pdflatex suppl.tex && pdflatex suppl.tex
```
