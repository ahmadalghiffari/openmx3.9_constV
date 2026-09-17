# Sohib et al., JJAP 64, 11SP11 (2025): reproduction note

## Scope

OpenMX 3.9.9 calculations were performed for the two geometries explicitly
discussed in the paper, `d=3.70 Ang` and `d=10.05 Ang`, at the zero-field
reference and `Delta f=2.6 eV`.  The calculations use PBE,
`C6.0-s2p2d1`, 250 Ry, a 20 Ang ESM cell, a 2.470 Ang graphene lattice,
and a `1x100x100` k grid (the ESM direction is stored as the first lattice
axis).  SOC was omitted in this initial reproduction; its effect for carbon
is expected to be much smaller than the numerical discrepancy discussed
below.

The released OpenMX parser reads `ESM.potential.diff`, not the longer
`ESM.potential.difference` seen in some example inputs.  The calculations
use `ESM.switch on4` (metal-vacuum-metal plus uniform field).

## Results at Delta f = 2.6 eV

| d (Ang) | transferred Q (electron/cell) | OpenMX Delta(VH/e), R-L (V) | physical Delta phi, R-L (V) | C/A (uF/cm2) |
|---:|---:|---:|---:|---:|
| 3.70 | 0.001885515 | -0.2212325 | +0.2212325 | 2.5844 |
| 10.05 | 0.005233962 | -0.5842369 | +0.5842369 | 2.7166 |

The primitive-cell area used in the conversion is 5.28353 Ang2.  The
transferred charge is the loss of electrons from the left graphene relative
to the zero-field calculation; the right graphene gains the same amount.

The paper explicitly quotes approximately 2.25 uF/cm2 for `d=10.05 Ang`
at `|Delta V|=0.6 V`.  The reproduced value at 0.584 V is 2.72 uF/cm2,
about 21% higher, but has the same scale and ordering.  Exact reproduction
will require the supplementary atomic coordinates/stacking and the authors'
ESM patch; these are not contained in the supplied PDF.

## Hartree-potential sign

OpenMX `dVHart` is the Hartree potential energy appearing in the one-electron
Hamiltonian.  For an electron of charge `-e`, the physical electrostatic
potential obeys

```text
phi = -V_H/e.
```

Consequently the OpenMX plane-averaged Hartree curve and the physical voltage
curve have opposite slopes.  In both reproduced cases OpenMX `V_H` decreases
from left to right, while physical `phi` increases from left to right.
Changing this sign does not change a capacitance reported as
`C=|Delta Q/Delta phi|`; it only changes a signed ratio or the stated voltage
polarity.  Thus the anticipated erratum can correct the potential/voltage
sign without invalidating the positive capacitance magnitudes.

## Files

The four `.dat` inputs, complete `.out` and `.log` files, Hartree-potential
cubes, and other OpenMX outputs are in this directory.  The plane-average
analysis uses `../capacitance_from_vhart.py`.
