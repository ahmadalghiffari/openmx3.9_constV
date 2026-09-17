#!/usr/bin/env python3
"""Compare ESM-on4 induced-Hartree diagnostics for 20 and 30 Ang cells."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt

from analyze_sohib_on4_components import (
    BOHR_TO_M, E_CHARGE, EPS0, HARTREE_TO_VOLT, cross, input_vesm,
    mulliken_populations, norm, periodic_interpolate, plane_average, read_cube,
)

DISTANCES = (10.05, 12.17, 14.29)
POSITIONS = {
    10.05: (0.332500000, 0.667500000),
    12.17: (0.297166667, 0.702833333),
    14.29: (0.261833333, 0.738166667),
}
TAGS = ("06", "12f", "18", "24", "30", "39")


def collect(distance):
    reference = mulliken_populations(f"on4_L30_d{distance:.2f}_vesm00.out")
    left, right = POSITIONS[distance]
    rows = []
    for tag in TAGS:
        stem = Path(f"on4_L30F_d{distance:.2f}_vesm{tag}")
        v_esm = input_vesm(str(stem) + ".dat")
        populations = mulliken_populations(str(stem) + ".out")
        delta = populations - reference
        q = 0.5 * (-delta[:2].sum() + delta[2:].sum())
        _, shape, steps, values = read_cube(str(stem) + ".vhart.cube")
        profile = plane_average(shape, values, 0)
        dvh_raw_signed = (
            periodic_interpolate(profile, right)
            - periodic_interpolate(profile, left)
        ) * HARTREE_TO_VOLT
        lattice = [tuple(shape[i] * x for x in steps[i]) for i in range(3)]
        area = norm(cross(lattice[1], lattice[2])) * BOHR_TO_M**2
        factor = E_CHARGE / area * 100.0
        # Poisson_ESM.c adds -0.5*V_ESM*(z-z1)/z1, whose slope is
        # -V_ESM/L.  Thus the signed right-minus-left external difference
        # is negative for positive V_ESM.
        dv_external = -v_esm * distance / 30.0
        dv_induced = abs(dvh_raw_signed - dv_external)
        rows.append({
            "distance_A": distance,
            "cell_A": 30.0,
            "edge_clearance_A": (30.0-distance)/2.0,
            "V_ESM_V": v_esm,
            "Q_e": q,
            "DeltaV_raw_signed_V": dvh_raw_signed,
            "DeltaV_external_layers_V": dv_external,
            "DeltaV_induced_V": dv_induced,
            "CH_induced_uF_cm2": factor*abs(q)/dv_induced,
            "Cparallel_eps0_over_d_uF_cm2": EPS0/(distance*1.0e-10)*100.0,
        })
    return rows


def main():
    rows = [row for d in DISTANCES for row in collect(d)]
    output = "sohib_L30_boundary_check.csv"
    with open(output, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    with open("sohib_stock_openmx_capacitance_components.csv") as stream:
        old = list(csv.DictReader(stream))
    old26 = {}
    for r in old:
        if float(r["V_ESM_V"]) != 2.6 or float(r["distance_A"]) not in DISTANCES:
            continue
        d = float(r["distance_A"])
        raw_abs = float(r["DeltaV_raw_Hartree_V"])
        # The stored old diagnostic used the wrong sign. Recover |Q|e/A from
        # Ctotal*|DeltaVraw| and subtract the signed external ramp correctly.
        charge_per_area = float(r["Ctotal_Sohib_uF_cm2"])*raw_abs
        induced = abs(-raw_abs - (-2.6*d/20.0))
        old26[d] = charge_per_area/induced
    new26 = {d: max((r for r in rows if r["distance_A"] == d),
                    key=lambda r: r["V_ESM_V"]) for d in DISTANCES}

    fig, ax = plt.subplots(figsize=(6.5, 4.6))
    ax.plot(DISTANCES, [old26[d] for d in DISTANCES], "o--", label="L=20 A")
    ax.plot(DISTANCES, [new26[d]["CH_induced_uF_cm2"] for d in DISTANCES],
            "s-", label="L=30 A (same external field)")
    ax.plot(DISTANCES,
            [new26[d]["Cparallel_eps0_over_d_uF_cm2"] for d in DISTANCES],
            "k^-.", label="epsilon0/d")
    ax.set_xlabel("Interelectrode distance d (A)")
    ax.set_ylabel("Capacitance per area (uF/cm2)")
    ax.set_title("ESM boundary check at the same external field")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig("sohib_L30_boundary_check.png", dpi=220)
    fig.savefig("sohib_L30_boundary_check.pdf")

    for d in DISTANCES:
        r = new26[d]
        print(f"d={d:.2f}: edge={r['edge_clearance_A']:.3f}, "
              f"Q={r['Q_e']:.9f}, raw={r['DeltaV_raw_signed_V']:.7f}, "
              f"induced={r['DeltaV_induced_V']:.7f}, "
              f"CH20={old26[d]:.6f}, CH30={r['CH_induced_uF_cm2']:.6f}, "
              f"eps0/d={r['Cparallel_eps0_over_d_uF_cm2']:.6f}")


if __name__ == "__main__":
    main()
