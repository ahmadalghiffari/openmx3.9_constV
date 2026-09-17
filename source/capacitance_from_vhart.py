#!/usr/bin/env python3
"""Calculate C=Q/delta(V_H) from an OpenMX Hartree-potential cube.

The potential is plane-averaged on planes normal to one cube lattice axis.
The two electrode positions are supplied as fractional coordinates along that
axis (0 <= f < 1).  Cube values and vectors are in atomic units.
"""

import argparse
import math
import re
from pathlib import Path

HARTREE_TO_VOLT = 27.211386245988
E_CHARGE = 1.602176634e-19
BOHR_TO_M = 5.29177210903e-11


def read_cube(path):
    with open(path, "r", encoding="utf-8") as f:
        f.readline()
        f.readline()
        fields = f.readline().split()
        natoms = abs(int(fields[0]))
        origin = tuple(map(float, fields[1:4]))
        shape = []
        steps = []
        for _ in range(3):
            fields = f.readline().split()
            shape.append(abs(int(fields[0])))
            steps.append(tuple(map(float, fields[1:4])))
        for _ in range(natoms):
            f.readline()
        values = [float(x) for line in f for x in line.split()]
    expected = shape[0] * shape[1] * shape[2]
    if len(values) != expected:
        raise ValueError(f"cube contains {len(values)} values; expected {expected}")
    return origin, shape, steps, values


def plane_average(shape, values, axis):
    n1, n2, n3 = shape
    sums = [0.0] * shape[axis]
    counts = [0] * shape[axis]
    p = 0
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                index = (i, j, k)[axis]
                sums[index] += values[p]
                counts[index] += 1
                p += 1
    return [s / n for s, n in zip(sums, counts)]


def periodic_interpolate(values, fraction):
    n = len(values)
    x = (fraction % 1.0) * n
    i = int(math.floor(x)) % n
    t = x - math.floor(x)
    return (1.0 - t) * values[i] + t * values[(i + 1) % n]


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def norm(a):
    return math.sqrt(sum(x*x for x in a))


def charge_from_log(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    matches = re.findall(r"Transferred charge:\s*([-+0-9.eE]+)", text)
    if not matches:
        raise ValueError("no 'Transferred charge' value found in log")
    return float(matches[-1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cube", help="OpenMX *.vhart.cube file")
    parser.add_argument("--axis", type=int, choices=(1, 2, 3), required=True,
                        help="cube/lattice axis normal to the electrode planes")
    parser.add_argument("--left", type=float, required=True,
                        help="left electrode fractional coordinate along axis")
    parser.add_argument("--right", type=float, required=True,
                        help="right electrode fractional coordinate along axis")
    charge = parser.add_mutually_exclusive_group(required=True)
    charge.add_argument("--charge", type=float, help="transferred charge in electrons")
    charge.add_argument("--log", help="OpenMX output containing the final summary")
    args = parser.parse_args()

    _, shape, steps, values = read_cube(args.cube)
    axis = args.axis - 1
    profile = plane_average(shape, values, axis)
    vleft_ha = periodic_interpolate(profile, args.left)
    vright_ha = periodic_interpolate(profile, args.right)
    # OpenMX dVHart is an electron potential energy.  The electrostatic
    # voltage therefore has the opposite sign: phi = -V_H/e.
    delta_vh_over_e = (vright_ha - vleft_ha) * HARTREE_TO_VOLT
    delta_phi = -delta_vh_over_e
    q_e = args.charge if args.charge is not None else charge_from_log(args.log)
    if abs(delta_phi) < 1.0e-15:
        raise ZeroDivisionError("Hartree potential difference is zero")
    capacitance = abs(q_e * E_CHARGE / delta_phi)

    other = [i for i in range(3) if i != axis]
    lattice_vectors = [tuple(shape[i] * x for x in steps[i]) for i in range(3)]
    area_m2 = norm(cross(lattice_vectors[other[0]], lattice_vectors[other[1]])) * BOHR_TO_M**2

    print(f"V_H(left)       = {vleft_ha: .12e} Hartree/e")
    print(f"V_H(right)      = {vright_ha: .12e} Hartree/e")
    print(f"Delta(V_H/e) R-L= {delta_vh_over_e: .12e} V")
    print(f"Delta phi R-L   = {delta_phi: .12e} V  (phi=-V_H/e)")
    print(f"Transferred Q   = {q_e: .12e} electron")
    print(f"Capacitance     = {capacitance: .12e} F")
    print(f"Capacitance     = {capacitance*1.0e18: .12e} aF")
    print(f"Area            = {area_m2: .12e} m^2")
    print(f"C / area        = {capacitance/area_m2: .12e} F/m^2")
    print(f"C / area        = {capacitance/area_m2*100.0: .12e} uF/cm^2")


if __name__ == "__main__":
    main()
