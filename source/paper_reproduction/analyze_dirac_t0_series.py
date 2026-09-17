#!/usr/bin/env python3
"""T=0 capacitances extracted solely from the electrostatic series relation."""

import csv
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from capacitance_from_vhart import (  # noqa: E402
    BOHR_TO_M, E_CHARGE, HARTREE_TO_VOLT, cross, norm,
    periodic_interpolate, plane_average, read_cube,
)

CASES = {"3.70": (0.40750, 0.59250), "10.05": (0.24875, 0.75125)}
COLORS = {"3.70": "tab:red", "10.05": "tab:blue"}


def last_float(pattern, text):
    values = re.findall(pattern, text)
    return float(values[-1]) if values else None


def collect(distance, left, right):
    rows = []
    for output in sorted(Path(".").glob(f"t0_d{distance}_vc*.std")):
        text = output.read_text(errors="replace")
        if "The calculation was normally finished." not in text:
            continue
        vc = last_float(r"Constraint difference Vc:\s*([-+0-9.eE]+)", text)
        q = last_float(r"Transferred charge:\s*([-+0-9.eE]+)", text)
        cube = output.with_suffix(".vhart.cube")
        if vc is None or q is None or not cube.exists():
            continue
        _, shape, steps, values = read_cube(cube)
        profile = plane_average(shape, values, 0)
        # OpenMX stores electron potential energy; physical phi_H has minus sign.
        dphi_signed = -(periodic_interpolate(profile, right)
                        - periodic_interpolate(profile, left))*HARTREE_TO_VOLT
        dphi = abs(dphi_signed)
        lattice = [tuple(shape[i]*x for x in steps[i]) for i in range(3)]
        area = norm(cross(lattice[1], lattice[2]))*BOHR_TO_M**2
        factor = E_CHARGE/area*100.0  # e/V/m2 -> microF/cm2
        ctot = factor*q/vc if vc != 0.0 else np.nan
        ch = factor*q/dphi if abs(dphi) > 1.0e-14 else np.nan
        # 1/Ctot = 1/CH + 1/CQ, equivalently CQ=Q/(Vc-dphi_H).
        residual = vc-dphi
        cq = factor*q/residual if abs(residual) > 1.0e-14 else np.nan
        rows.append({
            "distance_A": float(distance), "Vc_V": vc, "Q_e": q,
            "Delta_phi_H_V": dphi, "Vc_minus_Delta_phi_H_V": residual,
            "Ctotal_secant_uF_cm2": ctot, "CH_secant_uF_cm2": ch,
            "CQ_series_secant_uF_cm2": cq,
        })
    return rows


def main():
    data = {d: collect(d, *positions) for d, positions in CASES.items()}
    if any(len(rows) != 8 for rows in data.values()):
        raise SystemExit({d: len(rows) for d, rows in data.items()})
    rows = [row for case in data.values() for row in case]
    with open("graphene_dirac_t0_series.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(6.4, 4.7))
    for distance, case in data.items():
        x = np.array([r["Vc_V"] for r in case])
        y = np.array([r["CQ_series_secant_uF_cm2"] for r in case])
        ax.plot(np.r_[-x[::-1], x], np.r_[y[::-1], y], "o-",
                color=COLORS[distance], label=fr"$d={distance}$ $\AA$")
    ax.set_xlabel(r"Fixed potential difference $V_c$ (V)")
    ax.set_ylabel(r"Series-extracted secant $C_Q/A$ ($\mu$F cm$^{-2}$)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig("graphene_dirac_t0_cq_series.png", dpi=220)
    fig.savefig("graphene_dirac_t0_cq_series.pdf")

    for distance, case in data.items():
        print(f"d={distance} A")
        for row in case:
            print(f"  Vc={row['Vc_V']:.3f} Q={row['Q_e']:.8f} "
                  f"dphiH={row['Delta_phi_H_V']:.8f} "
                  f"CQ={row['CQ_series_secant_uF_cm2']:.8f}")


if __name__ == "__main__":
    main()
