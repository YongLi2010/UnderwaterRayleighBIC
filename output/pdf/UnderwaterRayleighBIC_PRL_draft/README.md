# Theory-only arXiv manuscript draft

This directory contains a PRL-style theory manuscript for the underwater
acoustic single-Rayleigh BIC.  It is limited to the infinite-period modal
model: no water-tank measurements, finite-element claims, finite-aperture
corrections, or fabricated experimental values are included.

## Files

- `main.tex` — PRL article source using `revtex4-2`.
- `suppl.tex` — supplemental derivation, convergence tables, caveats, and a
  clearly marked proposed future water-tank protocol.
- `ref.bib` — primary literature and foundational BIC references shared by
  both TeX sources.
- `figures/` — five main-text vector PDFs:
  `fig1_structure_topology.pdf`, `fig2_pole_linewidth_q.pdf`,
  `fig3_scattering_fields.pdf`, `fig4_cancellation_tolerance.pdf`, and
  `fig5_routing_device.pdf`.  The supplement also includes
  `figS1_experiment_protocol.pdf` and `figS2_controls_tolerance.pdf`;
  those panels are explicitly planning material, not measurements.

If a figure is not available, the TeX sources insert a framed placeholder so
that a syntax check can still run.  The current filenames above are the names
used by the manuscript.

## Compile

From this directory, use a standard BibTeX workflow for each source:

```sh
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
pdflatex suppl.tex && bibtex suppl && pdflatex suppl.tex && pdflatex suppl.tex
```

Alternatively, `latexmk -pdf main.tex` performs the same dependency cycle.
The manuscripts are designed for standard arXiv LaTeX services.  The delivered
PDFs were successfully rebuilt with bundled Tectonic 0.17.0: the main text is
7 pages and the supplement is 5 pages.  A standard TeX Live installation is
still recommended for the final arXiv upload check.

The numerical inputs, strict homogeneous-operator diagnostics, pole trend,
reciprocal routing values, and limitations are documented in \`main.tex\` and
trace to the modal-matching implementation in `../Ni2019_MATLAB/`.
