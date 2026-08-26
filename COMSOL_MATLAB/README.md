# COMSOL models for the rounded 180-kHz sample

This directory contains two solved COMSOL 6.4 models generated through
LiveLink for MATLAB:

- `rounded_180k_results/RayleighBIC_rounded_180k_eigenfrequency.mph`
- `rounded_180k_results/RayleighBIC_rounded_180k_scattering.mph`

## Fabrication geometry

All dimensions below are the rounded, two-decimal values used directly in the
COMSOL geometry.

| parameter | value (mm) |
|---|---:|
| period, `a` | 7.42 |
| large-groove width, `w1` | 4.41 |
| large-groove depth, `d1` | 5.15 |
| small-groove width, `w2` | 1.00 |
| small-groove depth, `d2` | 1.32 |
| intra-cell gap, `g12` | 1.00 |
| inter-cell edge gap, `g21` | 1.01 |
| total plate thickness | 10.00 |

The plate material node is stainless steel (`rho = 7850 kg m^-3`,
`E = 200 GPa`, `nu = 0.30`). Water is represented by `c = 1500 m s^-1`
and `rho = 1000 kg m^-3`.

## Physical scope

The solved files use the sound-hard water--steel interface. This is the direct
COMSOL counterpart of the modal-matching model used to identify the Rayleigh
BIC and is the appropriate baseline for comparing both methods. The explicit
10-mm steel domain and stainless-steel material are retained in the files.

A finite-elasticity calculation requires an Acoustic--Structure Boundary plus
a physically justified mechanical condition on the lower plate face (free,
clamped, or fixture impedance). That choice is experiment-specific and should
not be introduced silently because it adds Lamb-wave continua and can shift or
broaden the acoustic state.

## Solved results

- Selected eigenmode: `180.0463249 + 0.0003576 i kHz`
- Estimated radiation `Q`: `2.52e5`
- Driven scan: `179.60--180.40 kHz` in `0.02 kHz` steps at `7.10 deg`
- Maximum sampled groove response: `180.040 kHz`

The rounded geometry therefore realizes a near-BIC with a finite linewidth;
the strict machine-precision cancellation belongs to the unrounded optimized
geometry and is not claimed for this fabrication-rounded model.

## Reproduction

Start a COMSOL Multiphysics Server, connect LiveLink for MATLAB, and run:

```matlab
run('COMSOL_MATLAB/build_rounded_180k_models.m')
```

The field figure is generated from the exported complex COMSOL fields with:

```bash
python3 COMSOL_MATLAB/plot_comsol_fields.py
```

For display only, the field color limits use a symmetric 99.5-percentile clip
to prevent isolated sharp-corner FEM spikes from suppressing the spatial
pattern. The raw real and imaginary pressure arrays remain unchanged in CSV.
