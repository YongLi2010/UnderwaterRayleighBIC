# PRL submission evidence and figure plan
This plan is based on the current manuscript, the reproduction audit, and the completed analysis files. It preserves the central claim as a theory result: a strict Rayleigh-threshold homogeneous null in the lossless infinite-period two-groove modal model. It does not treat the reciprocal route as an isolator, the off-Rayleigh reflector as BIC-enhanced, or ideal similarity scaling as experimental environmental robustness.
## Recommended five-main-figure story
### Figure 1 — Geometry and the single-Rayleigh channel topology
Message: the design places one diffraction order exactly at the opening while a second order remains finite-flux open, so the threshold amplitude must be canceled in addition to ordinary radiation.
Suggested panels:
(a) Unit cell with \(a,d_1,d_2,w_1,w_2,g\), coordinate convention, and positive/negative incidence.
(b) Rayleigh lines for \(n=-1,0,+1\) in \((\kappa,\Omega)\), with the selected point and its order-reversed partner.
(c) Channel table or compact topology diagram showing \(n=0\) open, \(n=-1\) threshold, and \(n=+1\) evanescent for positive incidence.
(d) Normalized design point \(\kappa=0.0838486965529006\), \(\Omega=0.9161513034470994\), \(\theta=\pm5.251216456^\circ\).
The current source is arxiv_theory_paper/figures/fig1_structure_topology.pdf; its caption is main.tex Fig. \ref{fig:structure}. Verify that every dimension and channel label uses the final full-precision design cache Ni2019_MATLAB/results/StrictRayleighBIC_200kHz_min1mm.mat.
### Figure 2 — Strict null, full-operator verification, pole, and truncation
Message: the selected threshold pressure and the finite-flux open amplitude vanish at the same final root, and the pole continuation reaches that endpoint on the same outgoing branch.
Suggested panels:
(a) Same-geometry final-root pole trajectory in \((\kappa,\operatorname{Re}\Omega_p)\), with the \(n=-1\) Rayleigh line.
(b) \(Q\) versus \(|\Delta\kappa|\) from advanced_analysis/final_root_pole/final_root_pole_track.csv, annotating the local slope \(-1.99590867\) and the nearest finite value \(Q=2.0261339296\times10^{12}\) at \(|\Delta\kappa|=10^{-5}\).
(c) Full-operator endpoint evidence: reduced SVD ratio \(1.2603643174\times10^{-15}\), full complex-operator ratio \(7.7151137924\times10^{-16}\), raw residual \(1.3765668\times10^{-13}\), and unconstrained open/threshold amplitude errors.
(d) Reoptimized-root sequence versus fixed-high-root cross-check from advanced_analysis/extended_convergence/extended_root_sequence.csv and fixed_high_root_crosscheck.csv; explain that reoptimization tests existence of a sequence, whereas fixed-geometry checks expose truncation dependence.
Use advanced_analysis/final_root_pole/fig_final_root_pole.pdf and advanced_analysis/extended_convergence/figS_extended_convergence.pdf as data/graphic sources after redrawing into one consistent main figure. The old intermediate pole cache must not be silently used as the final-root continuation.
### Figure 3 — Homogeneous eigenfield localization
Message: the strict state is localized above the grating because all finite-flux and threshold exterior coefficients vanish; comparison states show how the exterior field changes away from the strict point.
Suggested panels:
(a) Strict BIC \(|p(x,y)|\) with a common normalization and the groove outlines.
(b) Controlled near-BIC homogeneous/quasi-mode field, explicitly labeled with its complex-pole or real-axis SVD-surrogate status.
(c) Ordinary leaky resonance field at a stated complex pole or clearly labeled real-axis surrogate.
(d) A compact quantitative inset: exterior-to-groove energy proxy, \(|A_{-1}|/\|A\|\), and \(|A_0|/\|A\|\) for all three states.
All states must use the same grid, branch convention, color scale, and coefficient normalization. Do not use the current driven panel (d) of fig3_scattering_fields.pdf as homogeneous-eigenfield evidence; the current figure manifest explicitly calls it a cached driven pressure field. The intended generator is advanced_analysis/robustness_fields.m. At this audit, field_states.csv, robustness_fields.mat, fig_homogeneous_fields.pdf, fig_dof_comparison.pdf, fig_fabrication_tolerance.pdf, and robustness_fields_report.md exist; they still require final independent numerical/visual inspection and manuscript integration.
### Figure 4 — Two independent radiation constraints and fabrication sensitivity
Message: two cavity degrees of freedom close two complex radiation phasors within the searched rectangular family, while fixed-geometry perturbations reopen the channels.
Suggested panels:
(a) The two-groove cell and the two complex constraints \(A_0^{\hom}=0\), \(A_{-1}^{\hom}=0\).
(b,c) Separate phasor diagrams for \(n=0\) and \(n=-1\), showing groove-1, groove-2, and total source vectors; quote closure magnitudes from advanced_analysis/reproduction_raw/14_phasors.txt.
(d) Bounded one-cavity comparison: \(\sigma_{\min}/\sigma_{\max}=0.1245882\) at (N,K)=(221,44) versus \(1.2604\times10^{-15}\) for the two-groove root; state the search family and do not call this a universal no-go theorem.
(e) Tolerance distributions for independent one-parameter and joint five-parameter perturbations at \(\pm1\%,\pm2\%,\pm5\%\), with raw residual and unconstrained threshold amplitude error separated.
The tolerance panels should only be included after tolerance_summary.csv has a report, reproducible seed, fixed/retuned definition, condition-number accounting, and a vector figure. The strict threshold amplitude is exactly zero by construction and is not a meaningful tolerance metric.
### Figure 5 — Reciprocal driven response, carefully qualified
Message: incidence reversal exchanges the signed anomalous order while preserving reciprocity and energy conservation; it does not produce isolation.
Suggested panels:
(a) \(\eta_{-1}(+\theta)\) and \(\eta_{+1}(-\theta)\) at the finite-truncation route point \(\theta=\pm6^\circ\), \(f=200.216874212\) kHz.
(b) Power-normalized scattering matrices in the ordered bases \((0,-1)\) and \((0,+1)\).
(c) Reciprocity and unitarity residuals; quote the direct values from reproduction_raw/13_reciprocity_direct.txt.
(d) A separate off-Rayleigh \(\eta_{\mathrm{target}}=0.9996185039\) point, visibly labeled ordinary and not BIC-enhanced, or move this panel to the supplement.
The final fixed-geometry route values are \(\eta_{-1}=0.3060611263\), \(\eta_0=0.6939388737\), and \(\theta_{\mathrm{out},-1}=-80.336964^\circ\) for positive incidence, with the order-reversed partner for negative incidence. Do not use the current hard-coded router-convergence table: suppl.tex Table \ref{tab:s-strict-convergence} is contradicted by direct fixed-frequency reruns. If a convergence inset is retained, it must show a documented joint retuning with all frequencies and condition numbers.
## Recommended supplementary figures and tables
S1. Modal-matching derivation and notation: the Floquet branch convention, groove cosine basis, projection matrix, full pole-free block operator, driven reduction, scaling, and power normalization. Use suppl.tex Eqs. \ref{eq:s-block}--\ref{eq:s-residuals}; include the exact threshold-tolerance rule.
S2. High-truncation convergence: the full figS_extended_convergence.pdf, with separate captions for reoptimized roots and fixed-high-geometry cross-checks. Include the complete CSV tables and full endpoint check.
S3. Conditional dimensional scaling and sound-speed control: advanced_analysis/scaling_environment/fig_scaling.pdf and fig_environment.pdf, with the exact law \(f_{\rm BIC}=\Omega_0c/(s a_0)\), the retuned \(c=1400\)--\(1600\) m/s results, and the fixed-device control showing that detuning is not a BIC.
S4. Detailed fabrication tolerance: vector distributions plus the complete tolerance_samples.csv and tolerance_summary.csv. Report independent and joint ensembles separately, with \(N,K\), random seed, number of draws, geometry-validity rules, raw residuals, unconstrained threshold/open amplitude errors, and condition numbers.
S5. Pole-track audit: full final-root CSV, complex branch convention, mode-overlap criterion, endpoint residual, and the finite-model nature of the \(Q\) fit. The old (121,15) cache may be retained only as a separately labeled historical comparison.
S6. Cavity-degree-of-freedom ablation: the search bounds, objective, transverse-mode treatment, direct single-cavity rerun, two-groove phasor sums, and dof_comparison.csv.
S7. Off-Rayleigh ordinary reflection: the full nonmonotonic truncation sequence and frequency/angle sensitivity at a converged modal order. This material is optional if the paper is shortened around the strict BIC.
S8. Proposed water-tank protocol: the current figS1_experiment_protocol.pdf and figS2_controls_tolerance.pdf may remain as planning schematics only. They must be labeled “proposed” and excluded from claims of validation.
## Prioritized simulations and checks before submission
P0 — Correct the evidence graph and regenerate all affected figures. Replace the old intermediate-root pole plot with the final-root track; remove or recompute the invalid router table; update \(\eta_{-1}\) to 0.3060611 if the finite-truncation route is retained; document the threshold classifier tolerance; ensure every panel states its \(N,K\), geometry, and whether it is homogeneous or driven.
P1 — Complete and validate the homogeneous-field and fabrication-tolerance deliverables. The current field_states.csv, robustness_fields.mat, fig_homogeneous_fields.pdf, fig_fabrication_tolerance.pdf, fig_dof_comparison.pdf, and robustness_fields_report.md must be independently inspected for numerical consistency, common normalization, and vector/raster status; add condition-number accounting to the tolerance summary; and distinguish exact strict zeros from unconstrained full-operator errors.
P2 — Validate modal convergence beyond the current finite sequence. Repeat the final-root strict and full-operator checks at further \((N,K)\), report root parameter drift, residuals, condition numbers, and branch/mode overlap. Add an independent implementation or full-wave check for the final geometry if the claim is intended to be more than a discretized modal result.
P3 — Quantify finite-device and nonideal effects: finite arrays with increasing period count, thermoviscous/viscous and thermal loss, water-property dispersion over 200 kHz, finite backing and fluid-structure compliance, source aperture, and disorder/fabrication perturbations. State whether the strict BIC becomes a high-Q resonance and how the threshold signature is operationally defined.
P4 — Strengthen the cavity-DOF conclusion only if desired. Expand the one-cavity search across a declared shape/mode family or provide a formal rank/independent-constraint argument. Until then retain the bounded wording “within the searched rectangular family.”
P5 — Decide whether the driven route belongs in the main narrative. If retained, perform joint frequency/geometry retuning at every truncation and report the condition number and energy closure. If not, move the route and the ordinary off-Rayleigh 99.96% result to the supplement so that the main five figures stay focused on strict BIC certification.
P6 — Complete submission metadata and reproducibility: replace placeholder author/affiliation fields, verify every BibTeX record/DOI/arXiv version, freeze MATLAB version and solver settings, archive scripts and machine-readable data, and make the exact figure-generation commands part of the release.
## Claim boundaries for the final draft
Safe: “The lossless infinite-period two-groove modal model has a strict Rayleigh-threshold homogeneous null at the reported normalized point, with both the threshold pressure amplitude and the finite-flux open amplitude numerically canceled.”
Safe with qualification: “A second independent cavity phase is needed within the declared rectangular search family,” “the pole continuation shows a local finite-model near-quadratic Q trend,” and “incidence reversal exchanges reciprocal signed diffraction orders.”
Not supported yet: “single-cavity impossibility,” “one-way/isolating device,” “robust experimental BIC,” “temperature-insensitive fixed device,” “universal Q exponent,” or “99.96% BIC-enhanced efficiency.”
