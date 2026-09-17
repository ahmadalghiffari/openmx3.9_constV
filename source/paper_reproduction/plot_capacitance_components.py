#!/usr/bin/env python3
"""Plot geometric, total, and extracted quantum capacitance for d=10.05 A."""

import csv
import glob
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from capacitance_from_vhart import (  # noqa: E402
    BOHR_TO_M, E_CHARGE, HARTREE_TO_VOLT, cross, norm,
    periodic_interpolate, plane_average, read_cube,
)

LEFT, RIGHT, AXIS = 0.24875, 0.75125, 0


def last_float(pattern, text):
    values = re.findall(pattern, text)
    return float(values[-1]) if values else None


def collect():
    files = glob.glob("graphene_d10.05_vc*_on1.log")
    files += glob.glob("graphene_d10.05_vc*_on1.std")
    rows, seen = [], set()
    for output in sorted(files):
        stem = output.removesuffix(".log").removesuffix(".std")
        if stem in seen:
            continue
        text = Path(output).read_text(errors="replace")
        if "The calculation was normally finished." not in text:
            continue
        q = last_float(r"Transferred charge:\s*([-+0-9.eE]+)", text)
        vc = last_float(r"Constraint difference Vc:\s*([-+0-9.eE]+)", text)
        if vc is None:
            old_half_vc = last_float(r"Final Vc:\s*([-+0-9.eE]+)", text)
            vc = 2.0*old_half_vc if old_half_vc is not None else None
        cube = stem + ".vhart.cube"
        if q is None or vc is None or not Path(cube).exists() or abs(vc) < 1e-14:
            continue
        _, shape, steps, values = read_cube(cube)
        profile = plane_average(shape, values, AXIS)
        vh_l = periodic_interpolate(profile, LEFT)
        vh_r = periodic_interpolate(profile, RIGHT)
        dphi_h = -(vh_r-vh_l)*HARTREE_TO_VOLT
        # Reject the known early-exit Vc=+0.04 point whose cube remained the
        # unconstrained Hartree grid despite a finite Mulliken transfer.
        if abs(dphi_h) < 1e-8 and abs(q) > 1e-10:
            continue
        lattice = [tuple(shape[i]*x for x in steps[i]) for i in range(3)]
        area = norm(cross(lattice[1], lattice[2]))*BOHR_TO_M**2
        vtotal = vc
        cgeom = abs(q)*E_CHARGE/(abs(dphi_h)*area)
        ctotal = abs(q)*E_CHARGE/(abs(vtotal)*area)
        inv_cq = 1.0/ctotal - 1.0/cgeom
        cq = 1.0/inv_cq if inv_cq > 0.0 else float("nan")
        rows.append({
            "Vc_eV": vc, "Vtotal_V": vtotal, "Q_e": q,
            "Vhart_physical_V": dphi_h,
            "Cgeom_uF_cm2": 100.0*cgeom,
            "Ctotal_uF_cm2": 100.0*ctotal,
            "Cq_uF_cm2": 100.0*cq,
            "output": output,
        })
        seen.add(stem)
    return sorted(rows, key=lambda r: r["Vtotal_V"])


def main():
    rows = collect()
    if not rows:
        raise SystemExit("No valid completed fixed-Vc calculations found")
    with open("graphene_d10.05_capacitance_components.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    # Add the symmetry origin in the same paper-style convention C_total(0)=0.
    plot_rows = list(rows)
    plot_rows.append({"Vtotal_V": 0.0, "Ctotal_uF_cm2": 0.0,
                      "Cgeom_uF_cm2": sum(r["Cgeom_uF_cm2"] for r in rows)/len(rows),
                      "Cq_uF_cm2": 0.0, "Q_e": 0.0})
    plot_rows.sort(key=lambda r: r["Vtotal_V"])

    fig, ax = plt.subplots(figsize=(6.4, 4.7))
    x = [r["Vtotal_V"] for r in plot_rows]
    ax.plot(x, [r["Ctotal_uF_cm2"] for r in plot_rows], "o-", label=r"$C_{\rm total}/A$")
    ax.plot(x, [r["Cgeom_uF_cm2"] for r in plot_rows], "s--", label=r"$C_{\rm geom}/A$ (Hartree)")
    ax.plot(x, [r["Cq_uF_cm2"] for r in plot_rows], "^--", label=r"extracted $C_Q/A$")
    ax.set_xlabel(r"Fixed constraint-potential difference $V_c$ (V)")
    ax.set_ylabel(r"Capacitance per area ($\mu$F cm$^{-2}$)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig("graphene_d10.05_capacitance_components.png", dpi=220)
    fig.savefig("graphene_d10.05_capacitance_components.pdf")

    fig, ax = plt.subplots(figsize=(6.4, 4.7))
    ax.plot([r["Vtotal_V"] for r in plot_rows],
            [r["Q_e"] for r in plot_rows], "o-")
    ax.set_xlabel(r"Fixed constraint-potential difference $V_c$ (V)")
    ax.set_ylabel(r"Transferred charge $Q$ (electron/cell)")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig("graphene_d10.05_q_vs_2vc.png", dpi=220)
    fig.savefig("graphene_d10.05_q_vs_2vc.pdf")


if __name__ == "__main__":
    main()
