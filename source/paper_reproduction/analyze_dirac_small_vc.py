#!/usr/bin/env python3
"""Differential-capacitance analysis of the dense small-Vc sweeps."""

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
    found = re.findall(pattern, text)
    return float(found[-1]) if found else None


def collect(distance, left, right):
    rows = []
    for output in sorted(Path(".").glob(f"dirac_d{distance}_vc*.std")):
        text = output.read_text(errors="replace")
        if "The calculation was normally finished." not in text:
            continue
        vc = last_float(r"Constraint difference Vc:\s*([-+0-9.eE]+)", text)
        q = last_float(r"Transferred charge:\s*([-+0-9.eE]+)", text)
        stem = str(output).removesuffix(".std")
        cube = Path(stem + ".vhart.cube")
        if vc is None or q is None or not cube.exists():
            continue
        _, shape, steps, values = read_cube(cube)
        profile = plane_average(shape, values, 0)
        dphi = -(periodic_interpolate(profile, right)
                 - periodic_interpolate(profile, left))*HARTREE_TO_VOLT
        lattice = [tuple(shape[i]*x for x in steps[i]) for i in range(3)]
        area = norm(cross(lattice[1], lattice[2]))*BOHR_TO_M**2
        rows.append((vc, q, abs(dphi), area, str(output)))
    return sorted(rows)


def differential_rows(distance, rows):
    vp = np.array([r[0] for r in rows])
    qp = np.array([r[1] for r in rows])
    hp = np.array([r[2] for r in rows])
    area = rows[0][3]
    voltage = np.r_[-vp[::-1], 0.0, vp]
    charge = np.r_[-qp[::-1], 0.0, qp]
    hartree = np.r_[-hp[::-1], 0.0, hp]
    dq_dv = np.gradient(charge, voltage, edge_order=2)
    dh_dv = np.gradient(hartree, voltage, edge_order=2)
    factor = E_CHARGE/area*100.0  # e/V/m2 -> microF/cm2
    ctotal = factor*dq_dv
    chartree = factor*dq_dv/dh_dv
    inverse_cq = 1.0/ctotal - 1.0/chartree
    cq = np.where(inverse_cq > 0.0, 1.0/inverse_cq, np.nan)
    result = []
    for i in range(len(voltage)):
        result.append({
            "distance_A": float(distance), "Vc_V": voltage[i],
            "Q_e": charge[i], "Delta_phi_H_V": hartree[i],
            "dQ_dVc_e_per_V": dq_dv[i],
            "dphiH_dVc": dh_dv[i],
            "Ctotal_diff_uF_cm2": ctotal[i],
            "CH_diff_uF_cm2": chartree[i],
            "CQ_diff_uF_cm2": cq[i],
        })
    return result


def linear_fit(rows):
    selected = [r for r in rows if 0.06 <= r["Vc_V"] <= 0.16]
    x = np.array([r["Vc_V"] for r in selected])
    y = np.array([r["CQ_diff_uF_cm2"] for r in selected])
    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope*x + intercept
    ss_res = np.sum((y-predicted)**2)
    ss_tot = np.sum((y-y.mean())**2)
    r2 = 1.0-ss_res/ss_tot if ss_tot > 0.0 else np.nan
    return slope, intercept, r2


def main():
    data = {}
    for distance, positions in CASES.items():
        raw = collect(distance, *positions)
        if len(raw) != 8:
            raise SystemExit(f"Expected 8 completed points for d={distance}, found {len(raw)}")
        data[distance] = differential_rows(distance, raw)

    all_rows = [r for rows in data.values() for r in rows]
    with open("graphene_dirac_small_vc_differential.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=all_rows[0].keys())
        writer.writeheader()
        writer.writerows(all_rows)

    fig, ax = plt.subplots(figsize=(6.4, 4.7))
    for distance, rows in data.items():
        ax.plot([r["Vc_V"] for r in rows], [r["Q_e"] for r in rows], "o-",
                color=COLORS[distance], label=fr"$d={distance}$ $\AA$")
    ax.set_xlabel(r"Fixed potential difference $V_c$ (V)")
    ax.set_ylabel(r"Transferred charge $Q$ (electron/cell)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig("graphene_dirac_small_vc_q.png", dpi=220)
    fig.savefig("graphene_dirac_small_vc_q.pdf")

    fig, ax = plt.subplots(figsize=(6.4, 4.7))
    for distance, rows in data.items():
        ax.plot([r["Vc_V"] for r in rows],
                [r["Ctotal_diff_uF_cm2"] for r in rows], "o-",
                color=COLORS[distance], label=fr"$d={distance}$ $\AA$")
    ax.set_xlabel(r"Fixed potential difference $V_c$ (V)")
    ax.set_ylabel(r"Differential $C_{\rm total}/A$ ($\mu$F cm$^{-2}$)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig("graphene_dirac_small_vc_ctotal_diff.png", dpi=220)
    fig.savefig("graphene_dirac_small_vc_ctotal_diff.pdf")

    fig, ax = plt.subplots(figsize=(6.4, 4.7))
    for distance, rows in data.items():
        slope, intercept, r2 = linear_fit(rows)
        ax.plot([r["Vc_V"] for r in rows], [r["CQ_diff_uF_cm2"] for r in rows],
                "o-", color=COLORS[distance],
                label=fr"$d={distance}$ $\AA$; outer linear $R^2={r2:.3f}$")
        print(f"d={distance}: CQ(0)={rows[len(rows)//2]['CQ_diff_uF_cm2']:.8f} "
              f"uF/cm2, outer slope={slope:.8f}, intercept={intercept:.8f}, R2={r2:.8f}")
    ax.axvspan(-0.026, 0.026, color="gray", alpha=0.12,
               label=r"$|V_c|<k_BT/e$ at 300 K")
    ax.set_xlabel(r"Fixed potential difference $V_c$ (V)")
    ax.set_ylabel(r"Series-extracted differential $C_Q/A$ ($\mu$F cm$^{-2}$)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8.5)
    fig.tight_layout()
    fig.savefig("graphene_dirac_small_vc_cq_diff.png", dpi=220)
    fig.savefig("graphene_dirac_small_vc_cq_diff.pdf")


if __name__ == "__main__":
    main()
