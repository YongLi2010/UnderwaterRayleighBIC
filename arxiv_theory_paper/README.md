# PRL-format theory manuscript

This directory contains the audited theory draft for a strict single-Rayleigh
BIC in a two-groove underwater acoustic metagrating.  The retained claim is
limited to an infinite, rigid, lossless periodic model.  No measurement,
finite-array, elastic-solid, or thermoviscous result is presented as evidence.

## Submission sources

- \`main.tex\` — four-figure PRL narrative using \`revtex4-2\`.
- \`suppl.tex\` — modal operator, endpoint audit, convergence, pole tracking,
  field definitions, bounded one-/two-cavity comparison, scaling,
  environmental controls, fabrication statistics, driven-scattering audit,
  and a clearly labeled proposed experiment.
- \`ref.bib\` — shared audited bibliography.
- \`figures/\` — referenced vector PDF figures:
  - main: \`fig1_structure_topology.pdf\`,
    \`fig2_strict_bic_evidence.pdf\`, \`fig3_homogeneous_fields.pdf\`,
    \`fig4_dof_tolerance.pdf\`;
  - supplement: \`figS1_extended_convergence.pdf\`,
    \`figS2_scaling_environment.pdf\`, \`figS3_tolerance_statistics.pdf\`,
    \`figS4_driven_scattering.pdf\`, \`figS5_proposed_experiment.pdf\`.

Additional older figure PDFs are retained only as development history and are
not referenced by the revised manuscript or included in the clean package.

## Compile

From this directory, a standard BibTeX workflow is:

\`\`\`sh
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
pdflatex suppl.tex && bibtex suppl && pdflatex suppl.tex && pdflatex suppl.tex
\`\`\`

The delivered files were also rebuilt with bundled Tectonic 0.17.0.  The
current output is six pages for the Letter, with the four evidence figures
integrated across pp. 2--5 and references completed on p. 6, plus eight pages
for the Supplemental Material.

## Evidence boundary

The final endpoint is the directly rerun \`(N,K)=(313,39)\` root.  The pole plot
uses a same-geometry continuation, not the superseded intermediate-root
cache.  The fixed-frequency near-threshold routing efficiency is explicitly
shown as truncation sensitive and is not a main-text device claim.  The
single-cavity comparison is a bounded searched-family result, not a universal
no-go theorem.
