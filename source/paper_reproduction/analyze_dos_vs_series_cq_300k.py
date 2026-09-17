#!/usr/bin/env python3
"""Compare 300 K layer-PDOS and fixed-Vc series quantum capacitances."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

E_CHARGE = 1.602176634e-19
KB_EV = 8.617333262145e-5
TEMPERATURE = 300.0
AREA_M2 = (2.470000000 * 2.139082748) * 1.0e-20
CONVERSION = E_CHARGE / AREA_M2 * 100.0  # (e/eV/cell)/A -> microF/cm2
COLORS = {3.7: "tab:red", 10.05: "tab:blue"}
PDOS_BASE = {3.7: "dos300_d3", 10.05: "dos300_d10"}


def read_atom_pdos(path):
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]


def read_layers(distance):
    base = PDOS_BASE[distance]
    energy = None
    atom_dos = []
    for atom in range(1, 5):
        e, dos = read_atom_pdos(Path(f"{base}.PDOS.Gaussian.atom{atom}"))
        if energy is None:
            energy = e
        elif not np.allclose(energy, e, atol=2.0e-6):
            raise ValueError("inconsistent PDOS energy meshes")
        atom_dos.append(dos)
    return energy, atom_dos[0] + atom_dos[1], atom_dos[2] + atom_dos[3]


def thermal_capacitance(energy, dos, mu):
    x = (energy-mu)/(2.0*KB_EV*TEMPERATURE)
    kernel = 1.0/(4.0*KB_EV*TEMPERATURE*np.cosh(x)**2)
    return CONVERSION*np.trapz(dos*kernel, energy)


def main():
    with open("graphene_dirac_small_vc_differential.csv") as stream:
        source = list(csv.DictReader(stream))

    output = []
    for distance in (3.7, 10.05):
        energy, dos_l, dos_r = read_layers(distance)
        selected = [r for r in source if abs(float(r["distance_A"])-distance) < 1e-6]
        for row in selected:
            vc = float(row["Vc_V"])
            dphi = float(row["Delta_phi_H_V"])
            # The non-Hartree voltage is shared symmetrically by the layers.
            vq = vc-dphi
            c_l = thermal_capacitance(energy, dos_l, +0.5*vq)
            c_r = thermal_capacitance(energy, dos_r, -0.5*vq)
            c_transfer = c_l*c_r/(c_l+c_r) if c_l+c_r > 0.0 else 0.0
            output.append({
                "distance_A": distance, "Vc_V": vc,
                "Delta_phi_H_V": dphi, "VQ_V": vq,
                "CQ_series_diff_uF_cm2": float(row["CQ_diff_uF_cm2"]),
                "CQ_DOS_left_uF_cm2": c_l,
                "CQ_DOS_right_uF_cm2": c_r,
                "CQ_DOS_transfer_uF_cm2": c_transfer,
            })

    with open("graphene_cq_dos_vs_series_300k.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=output[0].keys())
        writer.writeheader()
        writer.writerows(output)

    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    for distance in (3.7, 10.05):
        rows = [r for r in output if r["distance_A"] == distance]
        x = [r["Vc_V"] for r in rows]
        ax.plot(x, [r["CQ_series_diff_uF_cm2"] for r in rows], "o-",
                color=COLORS[distance],
                label=fr"series, $d={distance}$ $\AA$")
        ax.plot(x, [r["CQ_DOS_transfer_uF_cm2"] for r in rows], "--",
                color=COLORS[distance], linewidth=2,
                label=fr"layer-PDOS, $d={distance}$ $\AA$")
    ax.set_xlabel(r"Fixed potential difference $V_c$ (V)")
    ax.set_ylabel(r"Differential $C_Q/A$ ($\mu$F cm$^{-2}$)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8.5, ncol=2)
    fig.tight_layout()
    fig.savefig("graphene_cq_dos_vs_series_300k.png", dpi=220)
    fig.savefig("graphene_cq_dos_vs_series_300k.pdf")

    for distance in (3.7, 10.05):
        rows = [r for r in output if r["distance_A"] == distance]
        center = min(rows, key=lambda r: abs(r["Vc_V"]))
        edge = max(rows, key=lambda r: r["Vc_V"])
        print(f"d={distance}: at Vc=0 series={center['CQ_series_diff_uF_cm2']:.6f}, "
              f"DOS={center['CQ_DOS_transfer_uF_cm2']:.6f}; "
              f"at Vc={edge['Vc_V']:.2f} series={edge['CQ_series_diff_uF_cm2']:.6f}, "
              f"DOS={edge['CQ_DOS_transfer_uF_cm2']:.6f}")


if __name__ == "__main__":
    main()
