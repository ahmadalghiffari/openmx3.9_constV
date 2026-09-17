#!/usr/bin/env python3
"""Collect d=10.05 Ang cDFT/ESM-on1 runs and plot V versus C_total/A."""

import csv
import glob
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from capacitance_from_vhart import (  # noqa: E402
    BOHR_TO_M, E_CHARGE, HARTREE_TO_VOLT, cross, norm,
    periodic_interpolate, plane_average, read_cube,
)

LEFT = 0.24875
RIGHT = 0.75125
AXIS = 0
PATTERNS = ["graphene_d10.05_vc*_on1.log", "graphene_d10.05_vc*_on1.std"]


def final_charge(text):
    values = re.findall(r"Transferred charge:\s*([-+0-9.eE]+)", text)
    return float(values[-1]) if values else None


def collect():
    rows = []
    logs = []
    for pattern in PATTERNS:
        logs.extend(glob.glob(pattern))
    for log_name in sorted(set(logs)):
        text = Path(log_name).read_text(errors="replace")
        if "The calculation was normally finished." not in text:
            continue
        q = final_charge(text)
        stem = log_name.removesuffix(".log").removesuffix(".std")
        cube_name = stem + ".vhart.cube"
        if q is None or not Path(cube_name).exists():
            continue
        _, shape, steps, values = read_cube(cube_name)
        profile = plane_average(shape, values, AXIS)
        vh_l = periodic_interpolate(profile, LEFT)
        vh_r = periodic_interpolate(profile, RIGHT)
        dvh_rl = (vh_r - vh_l) * HARTREE_TO_VOLT
        dphi_rl = -dvh_rl
        # A run that exits immediately after switching on the constraint can
        # contain the pre-constraint Hartree grid.  Such a zero-voltage,
        # finite-charge point is not a valid converged capacitor state.
        if abs(dphi_rl) < 1.0e-8 and abs(q) > 1.0e-10:
            continue
        lattice = [tuple(shape[i] * x for x in steps[i]) for i in range(3)]
        area = norm(cross(lattice[1], lattice[2])) * BOHR_TO_M**2
        c_area = abs(q * E_CHARGE / dphi_rl) / area
        rows.append({"q_e": q, "voltage_V": dphi_rl,
                     "ct_F_m2": c_area, "ct_uF_cm2": 100.0*c_area,
                     "log": log_name})
    rows.sort(key=lambda r: abs(r["q_e"]))
    return rows


def main():
    positive = collect()
    if not positive:
        raise SystemExit("No completed sweep calculations found")

    # The geometry and boundary conditions are left-right symmetric.  Add the
    # negative-voltage branch by exact mirror reflection of the calculated
    # positive branch, and add the paper-style origin C_total(0)=0.
    rows = []
    for r in reversed(positive):
        m = dict(r)
        m["q_e"] = -r["q_e"]
        m["voltage_V"] = -r["voltage_V"]
        m["log"] = "mirror:" + r["log"]
        rows.append(m)
    rows.append({"q_e": 0.0, "voltage_V": 0.0, "ct_F_m2": 0.0,
                 "ct_uF_cm2": 0.0, "log": "symmetry origin"})
    rows.extend(positive)
    rows.sort(key=lambda r: r["voltage_V"])

    with open("graphene_d10.05_v_ct.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(6.2, 4.5))
    ax.plot([r["voltage_V"] for r in rows],
            [r["ct_uF_cm2"] for r in rows], "o-", lw=1.5, ms=5,
            label="cDFT + ESM on1")
    ax.set_xlabel(r"Physical voltage $\Delta\phi$ (V)")
    ax.set_ylabel(r"$C_\mathrm{total}/A$ ($\mu$F cm$^{-2}$)")
    ax.axvline(0.0, color="0.75", lw=0.8)
    ax.axhline(0.0, color="0.75", lw=0.8)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig("graphene_d10.05_v_ct.png", dpi=220)
    fig.savefig("graphene_d10.05_v_ct.pdf")

    # Also retain the primary Q-V relation used to compute the capacitance.
    fig, ax = plt.subplots(figsize=(6.2, 4.5))
    ax.plot([r["voltage_V"] for r in rows],
            [r["q_e"] for r in rows], "o-", lw=1.5, ms=5)
    ax.set_xlabel(r"Physical voltage $\Delta\phi$ (V)")
    ax.set_ylabel(r"Transferred charge $Q$ (electron/cell)")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig("graphene_d10.05_q_v.png", dpi=220)
    fig.savefig("graphene_d10.05_q_v.pdf")


if __name__ == "__main__":
    main()
