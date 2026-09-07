#!/usr/bin/env python
"""Project 01 - Cantilever beam: FE verification against beam theory + element/mesh convergence.

Pipeline
    Gmsh (structured hex or unstructured tet)  ->  CalculiX static step  ->  parse .frd / .dat
    ->  compare with Euler-Bernoulli and Timoshenko closed forms  ->  plots + PASS/FAIL checks

Run:  python run.py        (a few minutes on a laptop; results/ is rewritten, work/ is scratch)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fealib import analytic as an, ccx, mesh as fm, plotting

HERE = Path(__file__).resolve().parent
WORK, OUT = HERE / "work", HERE / "results"

# ---------------------------------------------------------------- problem definition (SI: m, N, Pa)
L, B, H = 0.400, 0.020, 0.020          # length, width (y), depth (z)
E, NU, RHO = 200e9, 0.30, 7850.0       # structural steel
P = 500.0                              # tip shear load, applied in -z as a uniform traction on x = L
SEC = an.rect_section(B, H)
G, KAPPA = an.shear_modulus(E, NU), an.shear_correction_rect(NU)
U_EB = an.cantilever_tip_deflection_eb(P, L, E, SEC["I"])
U_TIMO = an.cantilever_tip_deflection_timoshenko(P, L, E, SEC["I"], G, SEC["A"], KAPPA)
X_CHECK = 2 * H                        # stress checkpoint: two depths from the clamp (outside the Saint-Venant zone)

HEX_TYPES = ("C3D8", "C3D8R", "C3D8I", "C3D20", "C3D20R")
CASES = [(t, n) for t in HEX_TYPES for n in (1, 2, 4)] + [("C3D20R", 6)]
CASES += [("C3D4", n) for n in (1, 2, 4)] + [("C3D10", n) for n in (1, 2, 4)]
QUADRATIC = {"C3D20", "C3D20R", "C3D10"}


def build_mesh(etype: str, n: int) -> fm.Mesh:
    """n = elements through the depth (hex) or H/n target edge length (tet)."""
    order = 2 if etype in QUADRATIC else 1
    with fm.session("beam", order=order):
        if etype.startswith("C3D8") or etype.startswith("C3D20"):
            m = fm.box_hex(L, B, H, 20 * n, n, n, order=order)      # cubic elements of edge H/n
            return m.with_element_type({"C3D20" if order == 2 else "C3D8": etype})
        return fm.box_tet(L, B, H, H / n, order=order)


def solve(etype: str, n: int):
    tag = f"{etype}_n{n}"
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    m = build_mesh(etype, n)
    m.write_inp(wd / "mesh.inp")
    loads = m.traction_loads("XMAX", [0.0, 0.0, -P / (B * H)])     # consistent nodal loads, resultant = P
    deck = (
        f"** Cantilever {tag}: clamped at x=0, uniform shear traction on x=L\n"
        "*INCLUDE, INPUT=mesh.inp\n"
        + ccx.material_block("STEEL", E=E, nu=NU, rho=RHO)
        + "*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL\n*STEP\n*STATIC\n"
        + ccx.fmt_boundary("XMIN", 1, 3)
        + ccx.fmt_cload(loads)
        + "*NODE FILE\nU, S\n*NODE PRINT, NSET=XMIN, TOTALS=ONLY\nRF\n*END STEP\n"
    )
    (wd / "beam.inp").write_text(deck)
    res = ccx.run_ccx(wd / "beam.inp")
    frd, dat = ccx.read_frd(res.frd), ccx.read_dat(res.dat)
    disp, stress = frd.get("DISP"), frd.get("STRESS")

    u_tip = -disp.at(m.nsets["XMAX"])[:, 2].mean()
    rf = next(t for t in dat.totals if t["label"].startswith("force"))["values"]

    # top-fibre bending stress profile: average SXX over the nodes of z = H at each x station
    top = m.select_nodes(lambda x, y, z: abs(z - H) < 1e-9)
    prof = (pd.DataFrame({"x": np.round([m.nodes[i][0] for i in top], 9), "sxx": stress.at(top)[:, 0]})
            .groupby("x", as_index=False).mean())
    k = int(np.argmin(np.abs(prof.x.values - X_CHECK)))
    s_fe, x_fe = prof.sxx.values[k], prof.x.values[k]
    s_th = float(an.cantilever_bending_stress(x_fe, P, L, SEC["c"], SEC["I"]))

    # centre-line deflection (mean uz of each x station)
    allx = np.round([m.nodes[i][0] for i in disp.node_ids], 9)
    centre = pd.DataFrame({"x": allx, "uz": disp.values[:, 2]}).groupby("x", as_index=False).mean()

    row = dict(etype=etype, n=n, elements=m.n_elements, nodes=m.n_nodes, dofs=m.n_dofs(3),
               u_tip_mm=u_tip * 1e3, u_over_timoshenko=u_tip / U_TIMO, u_over_euler=u_tip / U_EB,
               x_check_m=x_fe, sxx_fe_MPa=s_fe / 1e6, sxx_theory_MPa=s_th / 1e6,
               sxx_err_pct=100 * (s_fe - s_th) / s_th, reaction_z_N=rf[2],
               equilibrium_err=abs(rf[2] - P) / P, solve_s=res.seconds)
    return row, prof, centre


def main() -> int:
    plotting.setup()
    OUT.mkdir(exist_ok=True)
    rows, profiles, centres = [], {}, {}
    for etype, n in CASES:
        row, prof, centre = solve(etype, n)
        rows.append(row)
        profiles[(etype, n)], centres[(etype, n)] = prof, centre
        print(f"{etype:7s} n={n}  dofs={row['dofs']:7d}  u/u_Timo={row['u_over_timoshenko']:.4f}  "
              f"sxx err={row['sxx_err_pct']:+.2f}%  eq.err={row['equilibrium_err']:.1e}  {row['solve_s']:.1f}s")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "convergence.csv", index=False)

    # ---------------- figure 1: convergence of tip deflection
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.axhline(1.0, color=plotting.INK2, lw=1.2, ls="--", label="Timoshenko (Euler-Bernoulli is 0.2 % lower)")
    ymax = 1.42
    for i, etype in enumerate(HEX_TYPES + ("C3D4", "C3D10")):
        d = df[df.etype == etype].sort_values("dofs")
        ax.plot(d.dofs, d.u_over_timoshenko.clip(upper=ymax), marker="o",
                ls="-" if etype.startswith(("C3D8", "C3D20")) else "--", color=plotting.PALETTE[i], label=etype)
    for _, r_ in df[df.u_over_timoshenko > ymax].iterrows():
        ax.annotate(f"{r_.etype}, n={r_.n}: {r_.u_over_timoshenko:.1f} x (hourglassing, off scale)", (r_.dofs, ymax),
                    textcoords="offset points", xytext=(10, -16), fontsize=8.5, color=plotting.INK2,
                    arrowprops=dict(arrowstyle="-", color=plotting.INK2, lw=0.8))
    ax.set_xscale("log")
    ax.set_ylim(0.4, 1.5)
    ax.set_xlabel("degrees of freedom")
    ax.set_ylabel("tip deflection / Timoshenko")
    ax.set_title("Cantilever tip deflection: element type and mesh convergence")
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.17), fontsize=8.5)
    plotting.save(fig, OUT / "convergence.png")

    # ---------------- figure 2: bending stress along the top fibre (finest hex mesh)
    key = ("C3D20R", 6)
    prof = profiles[key]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    xs = np.linspace(0, L, 200)
    ax.plot(xs * 1e3, an.cantilever_bending_stress(xs, P, L, SEC["c"], SEC["I"]) / 1e6, color=plotting.INK2, lw=1.5,
            ls="--", label="beam theory  M c / I")
    ax.plot(prof.x * 1e3, prof.sxx / 1e6, color=plotting.PALETTE[0], label=f"CalculiX {key[0]}, {df[(df.etype==key[0])&(df.n==key[1])].elements.item()} elements")
    ax.axvspan(0, X_CHECK * 1e3, color=plotting.PALETTE[3], alpha=0.12, lw=0)
    ax.text(X_CHECK * 1e3 * 0.5, ax.get_ylim()[1] * 0.92, "clamp\nzone", ha="center", va="top", fontsize=8.5, color=plotting.INK2)
    ax.set_xlabel("x from clamp (mm)")
    ax.set_ylabel("axial stress on top fibre (MPa)")
    ax.set_title("Bending stress: FE vs beam theory")
    ax.legend(loc="upper right")
    plotting.save(fig, OUT / "stress_profile.png")

    # ---------------- figure 3: deflection curve
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    for (etype, n), col in ((("C3D20R", 6), plotting.PALETTE[0]), (("C3D8", 1), plotting.PALETTE[1])):
        c = centres[(etype, n)]
        ax.plot(c.x * 1e3, c.uz * 1e3, color=col, lw=2.6, label=f"{etype}, n={n}")
    ax.plot(xs * 1e3, -an.cantilever_deflection_curve(xs, P, L, E, SEC["I"]) * 1e3, color=plotting.INK, lw=1.3, ls="--", label="Euler-Bernoulli", zorder=5)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("deflection (mm)")
    ax.set_title("Deflection curve: converged mesh vs a shear-locked coarse mesh")
    ax.legend(loc="lower left")
    plotting.save(fig, OUT / "deflection.png")

    # ---------------- objective checks
    best = df[(df.etype == "C3D20R") & (df.n == 6)].iloc[0]
    locked = df[(df.etype == "C3D8") & (df.n == 1)].iloc[0]
    checks = {
        "equilibrium: sum of clamp reactions equals applied load (all runs, < 1e-6)": bool((df.equilibrium_err < 1e-6).all()),
        "converged C3D20R tip deflection within 1 % of Timoshenko": bool(abs(best.u_over_timoshenko - 1) < 0.01),
        "converged C3D20R bending stress within 2 % of M c / I at x = 2H": bool(abs(best.sxx_err_pct) < 2),
        "quadratic elements converge monotonically (C3D20R)": bool(np.all(np.diff(df[df.etype == "C3D20R"].sort_values("dofs").u_over_timoshenko) >= -1e-4)),
        "shear locking reproduced: single-layer C3D8 under-predicts deflection by > 20 %": bool(locked.u_over_timoshenko < 0.8),
    }
    summary = {
        "inputs": dict(L=L, B=B, H=H, E=E, nu=NU, P=P, I=SEC["I"], kappa=KAPPA),
        "reference": dict(u_tip_euler_mm=U_EB * 1e3, u_tip_timoshenko_mm=U_TIMO * 1e3,
                          sigma_root_MPa=an.cantilever_bending_stress(0.0, P, L, SEC["c"], SEC["I"]) / 1e6),
        "best_run": {k: (float(v) if isinstance(v, (np.floating, float, int, np.integer)) else v) for k, v in best.items()},
        "checks": checks,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    (OUT / "convergence.md").write_text(df.to_markdown(index=False, floatfmt=".4g"))
    print("\nChecks:")
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
