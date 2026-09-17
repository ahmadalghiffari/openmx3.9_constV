#!/usr/bin/env python3
"""Reanalyse stock-OpenMX ESM-on4 data using Sohib's voltage definition."""

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

KB_EV = 8.617333262145e-5
T = 300.0
EPS0 = 8.8541878128e-12
CASES = {
    3.70: (0.40750, 0.59250), 5.82: (0.35450, 0.64550),
    7.93: (0.30175, 0.69825), 10.05: (0.24875, 0.75125),
    12.17: (0.19575, 0.80425), 14.29: (0.14275, 0.85725),
}
PDOS_BASE = {
    3.70: "dos300_d3", 5.82: "dos300_on4_d5",
    7.93: "dos300_on4_d7", 10.05: "dos300_d10",
    12.17: "dos300_on4_d12", 14.29: "dos300_on4_d14",
}


def mulliken_populations(path):
    text = Path(path).read_text(errors="replace")
    pattern = re.compile(
        r"^\s*([1-4])\s+C\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+([-+0-9.eE]+)\s+[-+0-9.eE]+\s*$",
        re.MULTILINE,
    )
    values = [(int(a), float(total)) for a, total in pattern.findall(text)]
    if len(values) < 4:
        raise ValueError(f"Mulliken table not found in {path}")
    return np.array([v for _, v in values[-4:]])


def layer_dos(distance):
    base = PDOS_BASE[distance]
    atoms = [np.loadtxt(f"{base}.PDOS.Gaussian.atom{i}") for i in range(1, 5)]
    energy = atoms[0][:, 0]
    return energy, atoms[0][:, 1]+atoms[1][:, 1], atoms[2][:, 1]+atoms[3][:, 1]


def dos_capacitance(energy, dos, mu, factor):
    x = (energy-mu)/(2.0*KB_EV*T)
    kernel = 1.0/(4.0*KB_EV*T*np.cosh(x)**2)
    return factor*np.trapz(dos*kernel, energy)


def input_vesm(path):
    text = Path(path).read_text(errors="replace")
    return float(re.search(r"ESM\.potential\.diff\s+([-+0-9.eE]+)", text).group(1))


def collect(distance, left, right):
    if distance in (3.70, 10.05):
        reference_stem = f"graphene_d{distance:.2f}_df0.0"
        stems = [Path(f"on4_d{distance:.2f}_vesm{tag}")
                 for tag in ("04", "08", "12", "16", "20")]
        stems.append(Path(f"graphene_d{distance:.2f}_df2.6"))
    else:
        reference_stem = f"on4_d{distance:.2f}_vesm00"
        stems = [Path(f"on4_d{distance:.2f}_vesm{tag}")
                 for tag in ("04", "08", "12", "16", "20", "26")]
    reference = mulliken_populations(reference_stem+".out")
    energy, dos_l, dos_r = layer_dos(distance)
    rows = []
    for stem in stems:
        dat = Path(str(stem)+".dat")
        out = Path(str(stem)+".out")
        cube = Path(str(stem)+".vhart.cube")
        v_esm = input_vesm(dat)
        populations = mulliken_populations(out)
        delta = populations-reference
        q = 0.5*(-delta[:2].sum()+delta[2:].sum())
        _, shape, steps, values = read_cube(cube)
        profile = plane_average(shape, values, 0)
        vl = periodic_interpolate(profile, left)
        vr = periodic_interpolate(profile, right)
        dvh_raw = (vr-vl)*HARTREE_TO_VOLT
        dv_raw = abs(dvh_raw)
        lattice = [tuple(shape[i]*x for x in steps[i]) for i in range(3)]
        area = norm(cross(lattice[1], lattice[2]))*BOHR_TO_M**2
        factor = E_CHARGE/area*100.0
        ctotal = factor*abs(q)/dv_raw

        # Sohib voltage is the full layer-to-layer raw Hartree difference.
        cl = dos_capacitance(energy, dos_l, +0.5*dv_raw, factor)
        cr = dos_capacitance(energy, dos_r, -0.5*dv_raw, factor)
        cq = cl*cr/(cl+cr)
        inverse_cgeom = 1.0/ctotal-1.0/cq
        cgeom_eff = 1.0/inverse_cgeom

        # Diagnostics not used in the Sohib series decomposition.
        # Poisson_ESM.c uses an external slope -V_ESM/L.
        dvext_layers = -v_esm*distance/20.0
        dvh_induced = dvh_raw-dvext_layers
        ch_induced = factor*abs(q)/abs(dvh_induced)
        cterminal = factor*abs(q)/abs(v_esm)
        cparallel = EPS0/(distance*1.0e-10)*100.0
        rows.append({
            "distance_A": distance, "V_ESM_V": v_esm,
            "Q_e": q, "DeltaV_raw_Hartree_V": dv_raw,
            "Ctotal_Sohib_uF_cm2": ctotal,
            "CQ_DOS_transfer_uF_cm2": cq,
            "Cgeom_effective_series_uF_cm2": cgeom_eff,
            "DeltaV_external_between_layers_V": dvext_layers,
            "DeltaV_induced_Hartree_V": abs(dvh_induced),
            "CH_induced_uF_cm2": ch_induced,
            "Cterminal_Q_over_VESM_uF_cm2": cterminal,
            "Cparallel_eps0_over_d_uF_cm2": cparallel,
        })
    return rows


def main():
    data = {d: collect(d, *positions) for d, positions in CASES.items()}
    rows = [r for case in data.values() for r in case]
    with open("sohib_stock_openmx_capacitance_components.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(2, 3, figsize=(12.2, 7.6), sharey=False)
    axes = axes.ravel()
    for ax, distance in zip(axes, CASES):
        case = data[distance]
        voltage = np.array([r["DeltaV_raw_Hartree_V"] for r in case])
        order = np.argsort(voltage)
        voltage = voltage[order]
        for key, marker, label in (
            ("Ctotal_Sohib_uF_cm2", "o-", r"$C_{\rm total}$ (Sohib)"),
            ("CQ_DOS_transfer_uF_cm2", "s--", r"$C_Q^{\rm DOS}$ (two layers)"),
            ("Cgeom_effective_series_uF_cm2", "^:", r"$C_{\rm geom}^{\rm eff}$ (series residual)"),
        ):
            y = np.array([r[key] for r in case])[order]
            ax.plot(voltage, y, marker, label=label)
        ax.set_title(fr"$d={distance:.2f}$ $\AA$")
        ax.set_xlabel(r"Sohib voltage $|\Delta V_H^{\rm raw}|$ (V)")
        ax.grid(alpha=0.25)
        ax.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel(r"Capacitance per area ($\mu$F cm$^{-2}$)")
    fig.tight_layout()
    fig.savefig("sohib_stock_openmx_capacitance_components.png", dpi=220)
    fig.savefig("sohib_stock_openmx_capacitance_components.pdf")

    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    distance = np.array(list(CASES))
    last = [max(data[d], key=lambda r: r["V_ESM_V"]) for d in distance]
    ax.plot(distance, [r["Ctotal_Sohib_uF_cm2"] for r in last], "o-",
            label=r"$C_{\rm total}$ (Sohib)")
    ax.plot(distance, [r["CQ_DOS_transfer_uF_cm2"] for r in last], "s--",
            label=r"$C_Q^{\rm DOS}$ (two layers)")
    ax.plot(distance, [r["Cgeom_effective_series_uF_cm2"] for r in last], "^:",
            label=r"$C_{\rm geom}^{\rm eff}$ (series residual)")
    ax.plot(distance, [r["CH_induced_uF_cm2"] for r in last], "d-.",
            label=r"$C_H^{\rm induced}$")
    ax.plot(distance, [r["Cparallel_eps0_over_d_uF_cm2"] for r in last], "k--",
            label=r"$\epsilon_0/d$")
    ax.set_xlabel(r"Interelectrode distance $d$ ($\AA$)")
    ax.set_ylabel(r"Capacitance per area ($\mu$F cm$^{-2}$)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig("sohib_stock_openmx_distance_dependence.png", dpi=220)
    fig.savefig("sohib_stock_openmx_distance_dependence.pdf")

    for distance, case in data.items():
        row = max(case, key=lambda r: r["V_ESM_V"])
        print(f"d={distance:.2f} A, VESM={row['V_ESM_V']:.1f}: "
              f"Vraw={row['DeltaV_raw_Hartree_V']:.7f}, Q={row['Q_e']:.9f}, "
              f"Ctotal={row['Ctotal_Sohib_uF_cm2']:.6f}, "
              f"CQdos={row['CQ_DOS_transfer_uF_cm2']:.6f}, "
              f"Cgeom_eff={row['Cgeom_effective_series_uF_cm2']:.6f}, "
              f"CHind={row['CH_induced_uF_cm2']:.6f}")


if __name__ == "__main__":
    main()
