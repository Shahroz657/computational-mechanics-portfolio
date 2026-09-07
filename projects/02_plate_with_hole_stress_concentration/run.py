#!/usr/bin/env python
"""Project 02 - Plate with a circular hole: stress-concentration factor sweep (2D plane stress).

Pipeline
    Gmsh quarter-plate with a distance-based size field around the hole  ->  CalculiX CPS8 plane
    stress  ->  Kt from nodal stresses  ->  compared with Peterson's chart fit (finite width) and
    Kirsch's infinite-plate field along the net section  ->  mesh convergence + PASS/FAIL checks.

Run:  python run.py     (about a minute; results/ is rewritten, work/ is scratch)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import gmsh
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fealib import analytic as an, ccx, mesh as fm, plotting

HERE = Path(__file__).resolve().parent
WORK, OUT = HERE / "work", HERE / "results"

# ---------------------------------------------------------------- problem definition (SI: m, N, Pa)
W, T = 0.100, 0.005            # plate width (x) and thickness
LH = 2.0 * W                   # modelled half-length (y): far-field boundary at 4 hole diameters or more
E, NU = 200e9, 0.30
S0 = 100e6                     # far-field tension in y
RATIOS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
HOLE_DIVS = 32                 # quadratic elements along the quarter arc (baseline)
CONV_RATIO, CONV_DIVS = 0.3, [8, 16, 32, 64]


def build(d_over_w: float, hole_divs: int = HOLE_DIVS):
    """Quarter plate [0, W/2] x [0, LH] minus a quarter hole of radius r at the origin."""
    r = 0.5 * d_over_w * W
    with fm.session("plate", order=2):
        occ = gmsh.model.occ
        rect = occ.addRectangle(0, 0, 0, W / 2, LH)
        disk = occ.addDisk(0, 0, 0, r, r)
        out, _ = occ.cut([(2, rect)], [(2, disk)])
        plate = out[0][1]
        occ.synchronize()
        tol = 1e-6 * W
        hole_tag = None
        for dim, tag in gmsh.model.getBoundary([(2, plate)], oriented=False):
            cx, cy, _ = occ.getCenterOfMass(dim, tag)
            if abs(cx) < tol:
                name = "SYM_X"
            elif abs(cy) < tol:
                name = "SYM_Y"
            elif abs(cy - LH) < tol:
                name = "TOP"
            elif abs(cx - W / 2) < tol:
                name = "RIGHT"
            else:
                name, hole_tag = "HOLE", tag
            gmsh.model.addPhysicalGroup(1, [tag], name=name)
        gmsh.model.addPhysicalGroup(2, [plate], name="PLATE")

        h_hole = (np.pi * r / 2) / hole_divs          # target edge length on the arc
        h_far = W / 12
        f = gmsh.model.mesh.field
        f.add("Distance", 1)
        f.setNumbers(1, "CurvesList", [hole_tag])
        f.setNumber(1, "Sampling", 400)
        f.add("Threshold", 2)
        f.setNumber(2, "InField", 1)
        f.setNumber(2, "SizeMin", h_hole)
        f.setNumber(2, "SizeMax", h_far)
        f.setNumber(2, "DistMin", 0.5 * r)
        f.setNumber(2, "DistMax", 5.0 * r)
        f.setAsBackgroundMesh(2)
        for opt in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
            gmsh.option.setNumber(opt, 0)
        gmsh.option.setNumber("Mesh.Algorithm", 8)            # Frontal-Delaunay for quads
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
        gmsh.model.mesh.generate(2)
        return fm.from_gmsh(2), r


def solve(d_over_w: float, hole_divs: int = HOLE_DIVS):
    tag = f"dW{d_over_w:.2f}_div{hole_divs}"
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    m, r = build(d_over_w, hole_divs)
    m.write_inp(wd / "mesh.inp")
    deck = (
        f"** Plate with hole, quarter model, d/W={d_over_w}, {hole_divs} elements on the arc\n"
        "*INCLUDE, INPUT=mesh.inp\n"
        + ccx.material_block("STEEL", E=E, nu=NU)
        + f"*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL\n{T}\n"
        "*STEP\n*STATIC\n"
        + ccx.fmt_boundary("SYM_X", 1, 1) + ccx.fmt_boundary("SYM_Y", 2, 2)
        + m.dload_lines("TOP", -S0)                                # negative pressure = tension
        + "*NODE FILE\nU, S\n*NODE PRINT, NSET=SYM_Y, TOTALS=ONLY\nRF\n*END STEP\n"
    )
    (wd / "plate.inp").write_text(deck)
    res = ccx.run_ccx(wd / "plate.inp")
    frd, dat = ccx.read_frd(res.frd), ccx.read_dat(res.dat)
    stress = frd.get("STRESS")
    syy = stress.column("SYY")
    hole = m.nsets["HOLE"]
    s_hole = stress.at(hole)[:, 1]
    k = int(np.argmax(s_hole))
    x_pk, y_pk, _ = m.nodes[hole[k]]
    s_max = float(s_hole[k])
    kt_gross = s_max / S0
    kt_net = kt_gross * (1.0 - d_over_w)
    rf = next(t for t in dat.totals if t["label"].startswith("force"))["values"]
    f_applied = S0 * (W / 2) * T
    # net-section profile along y = 0
    net = [n for n in m.nsets["SYM_Y"] if m.nodes[n][0] >= r - 1e-9]
    prof = pd.DataFrame({"x": [m.nodes[n][0] for n in net], "syy": stress.at(net)[:, 1]}).sort_values("x")
    row = dict(d_over_w=d_over_w, hole_divs=hole_divs, r_m=r, elements=m.n_elements, nodes=m.n_nodes,
               sigma_max_MPa=s_max / 1e6, kt_gross=kt_gross, kt_net=kt_net,
               kt_net_peterson=float(an.kt_plate_hole_net(d_over_w)),
               kt_err_pct=100 * (kt_net / float(an.kt_plate_hole_net(d_over_w)) - 1),
               peak_angle_deg=float(np.degrees(np.arctan2(y_pk, x_pk))),
               reaction_y_N=rf[1], equilibrium_err=abs(abs(rf[1]) - f_applied) / f_applied, solve_s=res.seconds)
    field = dict(zip(stress.node_ids.tolist(), (syy / S0).tolist()))
    return row, prof, (m, field)


def main() -> int:
    plotting.setup()
    OUT.mkdir(exist_ok=True)
    rows, profiles, fields = [], {}, {}
    for ratio in RATIOS:
        row, prof, fld = solve(ratio)
        rows.append(row)
        profiles[ratio], fields[ratio] = prof, fld
        print(f"d/W={ratio:.1f}  elems={row['elements']:6d}  Kt_net={row['kt_net']:.3f}  Peterson={row['kt_net_peterson']:.3f}  "
              f"err={row['kt_err_pct']:+.2f}%  peak@{row['peak_angle_deg']:.1f} deg  eq.err={row['equilibrium_err']:.1e}")
    sweep = pd.DataFrame(rows)
    sweep.to_csv(OUT / "kt_sweep.csv", index=False)

    conv_rows = []
    for divs in CONV_DIVS:
        row, _, _ = solve(CONV_RATIO, divs)
        conv_rows.append(row)
        print(f"convergence d/W={CONV_RATIO} divs={divs:3d} elems={row['elements']:6d}  Kt_net={row['kt_net']:.4f}")
    conv = pd.DataFrame(conv_rows)
    conv.to_csv(OUT / "kt_convergence.csv", index=False)

    # ---------------- figure 1: Kt vs d/W
    fig, ax = plt.subplots(figsize=(7.0, 4.3))
    rr = np.linspace(0.05, 0.65, 100)
    ax.plot(rr, an.kt_plate_hole_net(rr), color=plotting.INK2, lw=1.5, ls="--", label="Peterson chart fit (net section)")
    ax.plot(sweep.d_over_w, sweep.kt_net, "o", color=plotting.PALETTE[0], label="CalculiX CPS8, net-section Kt")
    ax.plot(rr, an.kt_plate_hole_gross(rr), color=plotting.INK2, lw=1.0, ls=":", label="Peterson, gross-section Kt")
    ax.plot(sweep.d_over_w, sweep.kt_gross, "s", color=plotting.PALETTE[1], label="CalculiX, gross-section Kt")
    for _, r_ in sweep.iterrows():
        ax.annotate(f"{r_.kt_net:.2f}", (r_.d_over_w, r_.kt_net), textcoords="offset points", xytext=(0, -14), ha="center", fontsize=8, color=plotting.INK2)
    ax.set_xlabel("hole diameter / plate width  d/W")
    ax.set_ylabel("stress concentration factor Kt")
    ax.set_title("Plate with a central hole in tension: Kt vs d/W")
    ax.legend(loc="upper left")
    plotting.save(fig, OUT / "kt_vs_dW.png")

    # ---------------- figure 2: net-section stress vs Kirsch (d/W = 0.1) and finite-width case (0.5)
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8), sharey=False)
    for ax, ratio in zip(axes, (0.1, 0.5)):
        prof = profiles[ratio]
        r = 0.5 * ratio * W
        xs = np.linspace(r, W / 2, 200)
        ax.plot(xs / r, an.kirsch_sigma_yy_along_x(xs, r, S0) / S0, color=plotting.INK2, lw=1.5, ls="--", label="Kirsch (infinite plate)")
        ax.plot(prof.x / r, prof.syy / S0, color=plotting.PALETTE[0], label="CalculiX (finite plate)")
        ax.set_xlabel("x / r along the net section")
        ax.set_ylabel("sigma_yy / sigma_0")
        ax.set_title(f"d/W = {ratio}")
        ax.legend(loc="upper right")
    fig.suptitle("Stress along the net section: infinite-plate theory vs finite-width FE", fontweight="bold")
    plotting.save(fig, OUT / "net_section_profile.png")

    # ---------------- figure 3: mesh convergence
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.plot(conv.elements, conv.kt_net, marker="o", color=plotting.PALETTE[0])
    for _, r_ in conv.iterrows():
        ax.annotate(f"{int(r_.hole_divs)} on arc", (r_.elements, r_.kt_net), textcoords="offset points", xytext=(6, -12), fontsize=8, color=plotting.INK2)
    ax.axhline(an.kt_plate_hole_net(CONV_RATIO), color=plotting.INK2, lw=1, ls="--")
    ax.text(conv.elements.max(), an.kt_plate_hole_net(CONV_RATIO), "  Peterson", va="center", fontsize=8.5, color=plotting.INK2)
    ax.set_xscale("log")
    ax.set_xlabel("elements")
    ax.set_ylabel("net-section Kt")
    ax.set_title(f"Mesh convergence at d/W = {CONV_RATIO}")
    plotting.save(fig, OUT / "kt_convergence.png")

    # ---------------- figure 4: stress field near the hole (d/W = 0.3)
    m, field = fields[0.3]
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    plotting.field_2d(ax, m, field, levels=30, mesh_lines=True, label="sigma_yy / sigma_0")
    ax.set_xlim(0, W / 2)
    ax.set_ylim(0, W / 2)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("sigma_yy / sigma_0 near the hole, d/W = 0.3 (quarter model, CPS8)")
    plotting.save(fig, OUT / "field_dW0.3.png")

    # ---------------- checks
    fine, mid = conv.iloc[-1], conv.iloc[-2]
    checks = {
        "equilibrium: symmetry-plane reaction equals applied edge load (all runs, < 1e-6)": bool((sweep.equilibrium_err < 1e-6).all()),
        "net-section Kt within 3 % of Peterson's chart fit for every d/W": bool((sweep.kt_err_pct.abs() < 3).all()),
        "peak stress sits on the net section (theta < 3 deg from the x axis)": bool((sweep.peak_angle_deg.abs() < 3).all()),
        "mesh convergence: Kt changes < 0.5 % between the two finest meshes": bool(abs(fine.kt_net / mid.kt_net - 1) < 0.005),
        "d/W -> 0 limit: gross Kt at d/W = 0.1 within 2 % of Kirsch's 3.0 corrected for width": bool(abs(sweep.iloc[0].kt_gross / an.kt_plate_hole_gross(0.1) - 1) < 0.02),
    }
    summary = {
        "inputs": dict(W=W, T=T, half_length=LH, E=E, nu=NU, sigma0=S0, hole_divs=HOLE_DIVS),
        "sweep": sweep.to_dict(orient="records"),
        "convergence": conv.to_dict(orient="records"),
        "checks": checks,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    (OUT / "kt_sweep.md").write_text(sweep[["d_over_w", "elements", "kt_gross", "kt_net", "kt_net_peterson", "kt_err_pct", "peak_angle_deg", "equilibrium_err"]].to_markdown(index=False, floatfmt=".4g"))
    print("\nChecks:")
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
