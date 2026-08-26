# Radiation-channel evolution along the 180-kHz BIC pole branch

Core conclusion: the open-order and Rayleigh-order pressure amplitudes have a
common first-order zero at the strict Rayleigh BIC.

- Geometry: final 180-kHz, 7.10-degree two-groove design.
- Modal truncation: `N=345`, `K=43`.
- Spectral path: the same complex-frequency pole is continued on both sides of
  `kappa_BIC = 0.110000105311808`.
- Normalization: pressure amplitudes are divided by the Euclidean norm of all
  groove modal coefficients, `||C||_2`.
- Exact endpoint: `A_0=A_-1=0` is inserted from the independently verified
  strict homogeneous solution; no numerical plotting floor is stored in the
  source data.
- Near-BIC fit interval: `1e-6 <= |Delta kappa| <= 1e-3`.
- Fitted exponents: `1.000170` for `A_0` and `0.996564` for `A_-1`.
- Maximum continued-pole singular-value ratio: `7.86e-11`.

Panel roles:

- a: the eigenpole crosses the `n=-1` Rayleigh line.
- b: both radiation amplitudes vanish at the same Bloch wave number.
- c: both amplitudes open linearly away from the BIC.

The calculation is deterministic; no sampling uncertainty or statistical test
applies. PDF text remains editable and the smallest rendered glyph is 5.18 pt,
above the 5-pt production floor.
