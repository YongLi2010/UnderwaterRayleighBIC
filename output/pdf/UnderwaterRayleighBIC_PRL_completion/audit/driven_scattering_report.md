# Driven-scattering audit at fixed final geometry

This audit independently reruns ni2019_modal_solver at the final two-groove geometry stored in Ni2019_MATLAB/results/StrictRayleighBIC_200kHz_min1mm.mat. No cached efficiency values or gallery tables are used. The geometry is held fixed while the modal truncation is varied.

The two operating points are intentionally separated:

- **Rayleigh-route point:** \(f=200216.874212\ {\rm Hz}\), \(\theta_i=+6^\circ\) and its signed partner \(-6^\circ\). The target order is \(n=-1\) for \(+6^\circ\) and \(n=+1\) for \(-6^\circ\).
- **Off-Rayleigh point:** \(f=202430\ {\rm Hz}\), \(\theta_i=+32.6328^\circ\), with \(n=-1\) as the target. This is an ordinary high-efficiency reflection state, not a BIC-enhanced state.

The complete solver-traced rows are in driven_scattering/driven_scattering_data.csv, and the MATLAB records and metadata are in driven_scattering/driven_scattering_audit.mat.

## Fixed-geometry Rayleigh route

The route response is strongly truncation-sensitive. At the specified frequency and geometry, \(\eta_{-1}\) is \(9.1843\times10^{-5}\), \(4.0379\times10^{-4}\), \(6.5929\times10^{-3}\), \(0.3060611\), and \(0.0571528\) for \((N,K)=(81,11),(121,15),(201,25),(313,39),(401,50)\), respectively. Thus the \(0.3060611\) value at \((313,39)\) is a finite-truncation illustrative operating point, not a converged fixed-geometry routing efficiency.

| \(N/K\) | \(\eta_{\rm target}\) | \(\eta_0\) | \(|1-\sum\eta|\) | condition number | \(\theta_{\rm out}\) (deg) |
|---:|---:|---:|---:|---:|---:|
| 81/11 | 0.0000918431 | 0.9999082 | \(1.11\times10^{-16}\) | \(1.9535\times10^4\) | -80.336964 |
| 121/15 | 0.0004037851 | 0.9995962 | \(9.99\times10^{-16}\) | \(4.6802\times10^4\) | -80.336964 |
| 201/25 | 0.0065929336 | 0.9934071 | \(1.33\times10^{-15}\) | \(2.6096\times10^5\) | -80.336964 |
| 313/39 | 0.3060611263 | 0.6939389 | \(3.73\times10^{-14}\) | \(2.5930\times10^6\) | -80.336964 |
| 401/50 | 0.0571527808 | 0.9428472 | \(3.16\times10^{-14}\) | \(1.4077\times10^6\) | -80.336964 |

The energy closure remains below \(9.9\times10^{-14}\) in the full audit, so the truncation sensitivity is not caused by a visible violation of the solver's lossless power balance. It is a response/convergence issue.

## Signed-angle reciprocity

At \((N,K)=(313,39)\), direct runs give

- \(+6^\circ\): \(\eta_{-1}=0.306061126315369\), \(\eta_0=0.693938873684669\), and \(\theta_{{\rm out},-1}=-80.336964032^\circ\);
- \(-6^\circ\): \(\eta_{+1}=0.306061126315479\), \(\eta_0=0.693938873684620\), and \(\theta_{{\rm out},+1}=+80.336964032^\circ\).

The signed-order partner differences are \(1.10\times10^{-13}\) in target efficiency, \(4.9\times10^{-14}\) in zero-order efficiency, and zero to printed precision after comparing opposite output signs. This is reciprocal signed-order exchange, not isolation or a one-way channel.

## Separate off-Rayleigh point

At \(f=202430\ {\rm Hz}\), \(\theta_i=32.6328^\circ\), the direct target efficiency increases through \((313,39)\) but decreases slightly at \((401,50)\): \(0.9779112\), \(0.9929755\), \(0.9986358\), \(0.9996289\), \(0.9999151\), and \(0.9996185\) over the six truncations. The final value is therefore reported as a high-efficiency ordinary state with a visible truncation-sensitivity caveat, not as a converged BIC-enhanced result.

The \((401,50)\) row has \(\eta_{-1}=0.999618503882140\), \(\eta_0=0.000381496117851\), \(|1-\sum\eta|=8.88\times10^{-15}\), condition number \(7.323\times10^3\), and \(\theta_{{\rm out},-1}=-32.626950^\circ\).

## Validation

Command:

    matlab -batch "run('advanced_analysis/driven_scattering_audit.m')"

The rerun completed without warnings or errors. The maximum energy closure error over all route, reciprocal, and off-Rayleigh rows was \(9.88\times10^{-14}\), below the \(10^{-8}\) acceptance threshold. The reciprocity partner at \((313,39)\) passed the \(10^{-8}\) target, zero-order, and signed-angle checks. The vector figure is one page, 501 x 375 pt, and pdfimages -list reports no raster images.

The four-panel figure uses a compact Helvetica vector layout. Panel (a) exposes the fixed-geometry route sensitivity; panel (b) shows reciprocal signed-order exchange at \(313/39\); panel (c) keeps the ordinary off-Rayleigh point separate; and panel (d) reports energy closure and conditioning alongside the route sweep.
