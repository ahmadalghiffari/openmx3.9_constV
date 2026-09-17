#!/usr/bin/env python3
"""Compare fixed-Vc cDFT capacitances for d=3.70 and 10.05 Angstrom."""

import csv
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from capacitance_from_vhart import (  # noqa: E402
    BOHR_TO_M, E_CHARGE, HARTREE_TO_VOLT, cross, norm,
    periodic_interpolate, plane_average, read_cube,
)

CASES = {
    "3.70": (0.40750, 0.59250),
    "10.05": (0.24875, 0.75125),
}
AXIS = 0


def last_float(pattern, text):
    values = re.findall(pattern, text)
    return float(values[-1]) if values else None


def collect(distance, left, right):
    rows, seen = [], set()
    outputs = list(Path(".").glob(f"graphene_d{distance}_vc0p*_on1.std"))
    outputs += list(Path(".").glob(f"graphene_d{distance}_vc0p*_on1.log"))
    for output in sorted(outputs):
        stem = str(output).removesuffix(".std").removesuffix(".log")
        if stem in seen:
            continue
        text = output.read_text(errors="replace")
        if "The calculation was normally finished." not in text:
            continue
        q = last_float(r"Transferred charge:\s*([-+0-9.eE]+)", text)
        vc = last_float(r"Constraint difference Vc:\s*([-+0-9.eE]+)", text)
        if vc is None:
            old_half_vc = last_float(r"Final Vc:\s*([-+0-9.eE]+)", text)
            vc = 2.0*old_half_vc if old_half_vc is not None else None
        cube = Path(stem + ".vhart.cube")
        if q is None or vc is None or vc <= 0.0 or not cube.exists():
            continue
        _, shape, steps, values = read_cube(cube)
        profile = plane_average(shape, values, AXIS)
        dphi = -(periodic_interpolate(profile, right)
                 - periodic_interpolate(profile, left))*HARTREE_TO_VOLT
        if abs(dphi) < 1.0e-10:
            continue
        lattice = [tuple(shape[i]*x for x in steps[i]) for i in range(3)]
        area = norm(cross(lattice[1], lattice[2]))*BOHR_TO_M**2
        vtotal = vc
        cgeom = abs(q)*E_CHARGE/(abs(dphi)*area)
        ctotal = abs(q)*E_CHARGE/(vtotal*area)
        inverse_cq = 1.0/ctotal - 1.0/cgeom
        cq = 1.0/inverse_cq if inverse_cq > 0.0 else float("nan")
        rows.append({
            "distance_A": float(distance), "Vc_eV": vc,
            "Vtotal_V": vtotal, "Q_e": q, "Delta_phi_H_V": dphi,
            "Cgeom_uF_cm2": 100.0*cgeom,
            "Ctotal_uF_cm2": 100.0*ctotal,
            "Cq_series_uF_cm2": 100.0*cq,
            "output": str(output),
        })
        seen.add(stem)
    return sorted(rows, key=lambda row: row["Vtotal_V"])


def symmetric(rows, field, origin):
    x = [-r["Vtotal_V"] for r in reversed(rows)] + [0.0]
    y = [r[field] for r in reversed(rows)] + [origin]
    x += [r["Vtotal_V"] for r in rows]
    y += [r[field] for r in rows]
    return x, y


def main():
    data = {d: collect(d, *positions) for d, positions in CASES.items()}
    if any(not rows for rows in data.values()):
        raise SystemExit("A completed sweep is missing")
    all_rows = [row for rows in data.values() for row in rows]
    with open("graphene_two_distance_cdft.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=all_rows[0].keys())
        writer.writeheader()
        writer.writerows(all_rows)

    colors = {"3.70": "tab:red", "10.05": "tab:blue"}

    fig, ax = plt.subplots(figsize=(6.4, 4.7))
    for d, rows in data.items():
        x, y = symmetric(rows, "Q_e", 0.0)
        ax.plot(x, y, "o-", color=colors[d], label=fr"$d={d}$ $\AA$")
    ax.set_xlabel(r"Fixed constraint-potential difference $V_c$ (V)")
    ax.set_ylabel(r"Transferred charge $Q$ (electron/cell)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig("graphene_two_distance_vc_q.png", dpi=220)
    fig.savefig("graphene_two_distance_vc_q.pdf")

    fig, ax = plt.subplots(figsize=(6.4, 4.7))
    for d, rows in data.items():
        x, y = symmetric(rows, "Ctotal_uF_cm2", 0.0)
        ax.plot(x, y, "o-", color=colors[d], label=fr"$d={d}$ $\AA$")
    ax.set_xlabel(r"Fixed constraint-potential difference $V_c$ (V)")
    ax.set_ylabel(r"Capacitance per area ($\mu$F cm$^{-2}$)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig("graphene_two_distance_fig4b_like.png", dpi=220)
    fig.savefig("graphene_two_distance_fig4b_like.pdf")

    fig, ax = plt.subplots(figsize=(6.4, 4.7))
    for d, rows in data.items():
        xg, yg = symmetric(rows, "Cgeom_uF_cm2",
                           sum(r["Cgeom_uF_cm2"] for r in rows)/len(rows))
        ax.plot(xg, yg, "s--", color=colors[d],
                label=fr"$d={d}$ $\AA$")
    ax.set_xlabel(r"Fixed constraint-potential difference $V_c$ (V)")
    ax.set_ylabel(r"Hartree $C_{\rm geom}/A$ ($\mu$F cm$^{-2}$)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig("graphene_two_distance_cgeom.png", dpi=220)
    fig.savefig("graphene_two_distance_cgeom.pdf")

    fig, ax = plt.subplots(figsize=(6.4, 4.7))
    for d, rows in data.items():
        x, y = symmetric(rows, "Cq_series_uF_cm2", 0.0)
        ax.plot(x, y, "o-", color=colors[d], label=fr"$d={d}$ $\AA$")
    ax.set_xlabel(r"Fixed constraint-potential difference $V_c$ (V)")
    ax.set_ylabel(r"Series-extracted $C_Q/A$ ($\mu$F cm$^{-2}$)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig("graphene_two_distance_fig4c_like_cq.png", dpi=220)
    fig.savefig("graphene_two_distance_fig4c_like_cq.pdf")


if __name__ == "__main__":
    main()
