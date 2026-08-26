# Fabrication parameters: 180-kHz Rayleigh BIC metagrating

The dimensions below are extracted directly from
`Ni2019_MATLAB/results/StrictRayleighBIC_180kHz_7p10deg_final.mat`.
The model is two-dimensional: the grooves are uniform along the out-of-plane
direction, open to water, and terminated by acoustically rigid walls.

| Quantity | Symbol | Nominal value |
|---|---:|---:|
| Lattice period | `a` | 7.416665789 mm |
| Large-groove width | `w1` | 4.408879373 mm |
| Large-groove depth | `d1` | 5.148350251 mm |
| Small-groove width | `w2` | 1.003212504 mm |
| Small-groove depth | `d2` | 1.324184849 mm |
| Large-to-small edge gap | `g12` | 1.002286956 mm |
| Small-to-next-large edge gap | `g21` | 1.002286956 mm |

The two edge-to-edge gaps are equal under periodic continuation.  With one
unit cell represented on `0 <= x < a`, the groove edges are:

| Feature | Left edge (mm) | Right edge (mm) |
|---|---:|---:|
| Large groove | 0.501143478 | 4.910022851 |
| Small groove | 5.912309807 | 6.915522311 |

The remaining 0.501143478 mm at each plotted unit-cell boundary joins across
the periodic boundary to form `g21`.

## Nominal operating point

- Water sound speed used in the design: 1500 m/s.
- Frequency: 180.000 kHz.
- Wavelength in water: 8.333333 mm.
- Bloch/incidence angle: 7.099662913 degrees.
- Normalized BIC coordinates: kappa = 0.110000105312 and Omega = 0.889999894688.

## Information required from a manufacturer

- Achievable tolerance and metrology for the 1.003 mm small groove and the
  1.002 mm gaps (request quotations for ±0.02 mm and ±0.05 mm).
- Minimum internal corner and bottom-corner radius.
- Surface roughness inside both grooves.
- Flatness and parallelism of the reflecting face.
- Candidate plate material and the residual back-wall thickness below the
  5.148 mm groove.

The theoretical model does not fix the finite number of periods, the
out-of-plane groove length, mounting frame, or total plate thickness.  Those
dimensions must be selected together with the source aperture and scan
distance.  The as-built geometry should be measured and inserted back into
the modal-matching model before setting the experimental frequency and angle.
