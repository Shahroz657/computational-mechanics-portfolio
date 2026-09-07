#!/usr/bin/env python
"""Project 05 - L-bracket design optimisation: CadQuery -> Gmsh -> CalculiX -> response surface -> SciPy.

Objective: minimum mass. Constraints: safety factor on fillet von Mises stress >= 2.0 against yield,
tip deflection <= DEFL_MAX. Design variables: plate thickness t and inner fillet radius R.
Method: a 7 x 5 design-of-experiments of full FE runs, quadratic response surfaces (in log space for
stress and deflection), constrained minimisation with SLSQP, then FE verification of the optimum and a
mesh-convergence check that separates the converged fillet stress from the singular stress at the
edge of the fixed face.

Unit system: mm - N - MPa (E = 200 000 MPa). Mass from the exact CAD volume.
Run:  python run.py     (about 5-8 minutes; results/ rewritten, work/ scratch)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import gmsh
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.optimize import minimize

from bracket import HV, LA, WB, export_step, export_svg, mass_kg
from fealib import ccx, mesh as fm, plotting

HERE = Path(__file__).resolve().parent
WORK, OUT = HERE / "work", HERE / "results"

# ---------------------------------------------------------------- requirements (mm - N - MPa)
E, NU, RHO = 200000.0, 0.30, 7850.0
S_Y = 355.0                       # S355 structural steel yield strength
F = 600.0                         # downward load on the arm end face
SF_REQ = 2.0                      # required safety factor on fillet von Mises stress
DEFL_MAX = 0.60                   # mm, allowable deflection of the loaded end
T_RANGE, R_RANGE = (5.0, 12.0), (2.0, 10.0)
DOE_T, DOE_R = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0], [2.0, 4.0, 6.0, 8.0, 10.0]
BASELINE = (10.0, 2.0)            # "thick plate, small fillet" reference design


def build_mesh(t: float, R: float, size: float = 6.0, fillet_div: int = 6) -> fm.Mesh:
    step = export_step(t, R, WORK / f"geom_t{t:g}_R{R:g}.step")
    with fm.session("bracket", order=2):
        occ = gmsh.model.occ
        occ.importShapes(str(step))
        occ.synchronize()
        tol = 1e-3
        fixed, load, fillet, skin = [], [], [], []
        for dim, tag in gmsh.model.getEntities(2):
            skin.append(tag)
            com = occ.getCenterOfMass(dim, tag)
            bb = gmsh.model.getBoundingBox(dim, tag)
            if abs(com[0]) < tol:
                fixed.append(tag)
            elif abs(com[0] - (t + LA)) < tol:
                load.append(tag)
            elif bb[0] > t - tol and bb[3] < t + R + tol and bb[2] > t - tol and bb[5] < t + R + tol:
                fillet.append(tag)
        assert len(fixed) == 1 and len(load) == 1 and len(fillet) == 1, (fixed, load, fillet)
        for name, tags in (("FIXED", fixed), ("LOAD", load), ("FILLET", fillet), ("SKIN", skin)):
            gmsh.model.addPhysicalGroup(2, tags, name=name)
        gmsh.model.addPhysicalGroup(3, [tag for _, tag in gmsh.model.getEntities(3)], name="BRACKET")
        f = gmsh.model.mesh.field
        f.add("Distance", 1)
        f.setNumbers(1, "SurfacesList", fillet)
        f.setNumber(1, "Sampling", 100)
        f.add("Threshold", 2)
        f.setNumber(2, "InField", 1)
        f.setNumber(2, "SizeMin", max(R / fillet_div, 0.5))
        f.setNumber(2, "SizeMax", size)
        f.setNumber(2, "DistMin", 0.5 * R)
        f.setNumber(2, "DistMax", 2.0 * R)
        f.setAsBackgroundMesh(2)
        gmsh.option.setNumber("Mesh.MeshSizeMax", size)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 16)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)          # HXT
        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.model.mesh.generate(3)
        return fm.from_gmsh(3)


def solve(t: float, R: float, size: float = 6.0, fillet_div: int = 6):
    tag = f"t{t:g}_R{R:g}_s{size:g}_d{fillet_div}"
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    m = build_mesh(t, R, size, fillet_div)
    m.write_inp(wd / "mesh.inp")
    loads = m.traction_loads("LOAD", [0.0, 0.0, -F / (WB * t)])
    deck = (
        f"** L-bracket {tag}: back face fixed, {F} N down on the arm end face (mm-N-MPa)\n"
        "*INCLUDE, INPUT=mesh.inp\n" + ccx.material_block("STEEL", E=E, nu=NU)
        + "*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL\n*STEP\n*STATIC\n"
        + ccx.fmt_boundary("FIXED", 1, 3) + ccx.fmt_cload(loads)
        + "*NODE FILE\nU, S\n*NODE PRINT, NSET=FIXED, TOTALS=ONLY\nRF\n*END STEP\n"
    )
    (wd / "bracket.inp").write_text(deck)
    res = ccx.run_ccx(wd / "bracket.inp")
    frd, dat = ccx.read_frd(res.frd), ccx.read_dat(res.dat)
    disp, stress = frd.get("DISP"), frd.get("STRESS")
    vm = ccx.von_mises(stress.values)
    vm_by = dict(zip(stress.node_ids.tolist(), vm.tolist()))
    s_fillet = max(vm_by[n] for n in m.nsets["FILLET"])
    s_fixed = max(vm_by[n] for n in m.nsets["FIXED"])
    defl = -disp.at(m.nsets["LOAD"])[:, 2].mean()
    rf = next(x for x in dat.totals if x["label"].startswith("force"))["values"]
    row = dict(t=t, R=R, size=size, fillet_div=fillet_div, elements=m.n_elements, nodes=m.n_nodes,
               mass_kg=mass_kg(t, R, RHO), sigma_fillet_MPa=s_fillet, sigma_fixed_edge_MPa=s_fixed,
               sigma_global_max_MPa=float(vm.max()), SF=S_Y / s_fillet, deflection_mm=defl,
               reaction_z_N=rf[2], equilibrium_err=abs(rf[2] - F) / F, solve_s=res.seconds)
    return row, (m, vm_by, disp)


# ---------------------------------------------------------------- response surfaces
def _features(t, R):
    t, R = np.asarray(t, float), np.asarray(R, float)
    return np.stack([np.ones_like(t), t, R, t * t, R * R, t * R], axis=-1)


def fit_quadratic(df: pd.DataFrame, col: str, log: bool):
    y = np.log(df[col].values) if log else df[col].values
    X = _features(df.t.values, df.R.values)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coef
    r2 = 1.0 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
    back = np.exp(pred) if log else pred
    return coef, float(r2), float(np.max(np.abs(back / df[col].values - 1.0)))


def predict(coef, t, R, log: bool):
    v = _features(t, R) @ coef
    return np.exp(v) if log else v


def stress_3d_figure(m: fm.Mesh, vm_by: dict, disp, path, scale: float = 40.0, title: str = ""):
    u = disp.as_dict()
    polys, vals = [], []
    for et, row in m.boundary["SKIN"]:
        corners = row[:3] if et in (2, 9) else row[:4]
        polys.append([np.array(m.nodes[n]) + scale * u[n] for n in corners])
        vals.append(np.mean([vm_by[n] for n in corners]))
    fig = plt.figure(figsize=(7.6, 5.4))
    ax = fig.add_subplot(111, projection="3d")
    coll = Poly3DCollection(polys, edgecolor="none")
    coll.set_array(np.array(vals))
    coll.set_cmap(plotting.sequential_cmap())
    coll.set_clim(0, S_Y / SF_REQ)
    ax.add_collection3d(coll)
    ax.set_xlim(-5, LA + 15)
    ax.set_ylim(-5, WB + 5)
    ax.set_zlim(-25, HV + 5)
    ax.set_box_aspect((LA + 20, WB + 10, HV + 30))
    ax.view_init(elev=26, azim=-52)
    ax.set_axis_off()
    cb = fig.colorbar(coll, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label(f"von Mises stress (MPa), capped at allowable {S_Y / SF_REQ:.0f} MPa")
    ax.set_title(title, fontsize=10)
    plotting.save(fig, path)


def main() -> int:
    plotting.setup()
    OUT.mkdir(exist_ok=True)
    WORK.mkdir(exist_ok=True)

    # ---------------- design of experiments
    rows = []
    for t in DOE_T:
        for R in DOE_R:
            row, _ = solve(t, R)
            rows.append(row)
            print(f"DOE t={t:4.1f} R={R:4.1f}  elems={row['elements']:6d}  fillet={row['sigma_fillet_MPa']:6.1f} MPa  "
                  f"SF={row['SF']:.2f}  defl={row['deflection_mm']:.3f} mm  mass={row['mass_kg']:.3f} kg  {row['solve_s']:.1f}s")
    doe = pd.DataFrame(rows)
    doe.to_csv(OUT / "doe.csv", index=False)

    c_s, r2_s, res_s = fit_quadratic(doe, "sigma_fillet_MPa", log=True)
    c_d, r2_d, res_d = fit_quadratic(doe, "deflection_mm", log=True)
    c_m, r2_m, res_m = fit_quadratic(doe, "mass_kg", log=False)
    print(f"RSM fit: stress R2={r2_s:.4f} (max resid {100*res_s:.1f} %), deflection R2={r2_d:.4f} (max resid {100*res_d:.1f} %), mass R2={r2_m:.5f}")

    # ---------------- constrained minimisation on the surrogate
    cons = [{"type": "ineq", "fun": lambda x: S_Y / predict(c_s, x[0], x[1], True) - SF_REQ},
            {"type": "ineq", "fun": lambda x: DEFL_MAX - predict(c_d, x[0], x[1], True)}]
    best = None
    for t0 in (6.0, 8.0, 10.0):
        for R0 in (3.0, 6.0, 9.0):
            r = minimize(lambda x: predict(c_m, x[0], x[1], False), [t0, R0], method="SLSQP",
                         bounds=[T_RANGE, R_RANGE], constraints=cons)
            feasible = all(c["fun"](r.x) > -1e-6 for c in cons)
            if r.success and feasible and (best is None or r.fun < best.fun):
                best = r
    t_opt, R_opt = float(best.x[0]), float(best.x[1])
    print(f"surrogate optimum: t={t_opt:.2f} mm, R={R_opt:.2f} mm, mass={best.fun:.4f} kg, "
          f"SF_pred={S_Y/predict(c_s, t_opt, R_opt, True):.3f}, defl_pred={predict(c_d, t_opt, R_opt, True):.3f} mm")

    # ---------------- FE verification of the optimum, with one corrective trim if needed
    ver, _ = solve(round(t_opt, 2), round(R_opt, 2))
    print(f"FE at optimum: SF={ver['SF']:.3f}, defl={ver['deflection_mm']:.3f} mm, mass={ver['mass_kg']:.4f} kg")
    trimmed = False
    if ver["SF"] < SF_REQ or ver["deflection_mm"] > DEFL_MAX:
        trimmed = True
        factor = max((SF_REQ / ver["SF"]) ** 0.5, (ver["deflection_mm"] / DEFL_MAX) ** (1 / 3), 1.0)   # sigma ~ 1/t^2, defl ~ 1/t^3
        t_fix = round(min(T_RANGE[1], ver["t"] * factor * 1.005), 2)
        ver, _ = solve(t_fix, round(R_opt, 2))
        print(f"trimmed design t={t_fix} mm: SF={ver['SF']:.3f}, defl={ver['deflection_mm']:.3f} mm, mass={ver['mass_kg']:.4f} kg")
    base, _ = solve(*BASELINE)
    print(f"baseline t={BASELINE[0]} R={BASELINE[1]}: SF={base['SF']:.3f}, defl={base['deflection_mm']:.3f} mm, mass={base['mass_kg']:.4f} kg")

    # ---------------- mesh convergence at the verified design: fillet stress vs singular fixed-edge stress
    conv_rows = []
    fields = None
    for size, div in ((8.0, 4), (6.0, 6), (4.0, 8), (3.0, 10)):
        row, fld = solve(ver["t"], ver["R"], size, div)
        conv_rows.append(row)
        fields = fld
        print(f"convergence size={size} div={div}: nodes={row['nodes']:7d} fillet={row['sigma_fillet_MPa']:.1f} MPa  fixed-edge={row['sigma_fixed_edge_MPa']:.1f} MPa")
    conv = pd.DataFrame(conv_rows)
    conv.to_csv(OUT / "convergence.csv", index=False)

    # ---------------- figures
    export_svg(ver["t"], ver["R"], OUT / "bracket_optimum.svg")
    stress_3d_figure(*fields, OUT / "stress_optimum.png", title=f"Optimised bracket t = {ver['t']:g} mm, R = {ver['R']:g} mm, deformation x40")

    tt, rr = np.meshgrid(np.linspace(*T_RANGE, 120), np.linspace(*R_RANGE, 120))
    sf = S_Y / predict(c_s, tt, rr, True)
    dd = predict(c_d, tt, rr, True)
    mm = predict(c_m, tt, rr, False)
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    cs = ax.contourf(tt, rr, mm, levels=18, cmap=plotting.sequential_cmap(), alpha=0.9)
    fig.colorbar(cs, ax=ax, label="mass (kg), response surface")
    ax.contour(tt, rr, sf, levels=[SF_REQ], colors=[plotting.PALETTE[7]], linewidths=2)
    ax.contour(tt, rr, dd, levels=[DEFL_MAX], colors=[plotting.PALETTE[1]], linewidths=2, linestyles="--")
    infeasible = (sf < SF_REQ) | (dd > DEFL_MAX)
    ax.contourf(tt, rr, infeasible.astype(float), levels=[0.5, 1.5], colors=["#ffffff"], alpha=0.55)
    ax.grid(False)
    ax.plot(doe.t, doe.R, "o", ms=4, color=plotting.INK2, label="FE design points (DOE)")
    ax.plot(BASELINE[0], BASELINE[1], "s", ms=9, color=plotting.INK, label=f"baseline {base['mass_kg']:.3f} kg, SF {base['SF']:.2f}")
    ax.plot(ver["t"], ver["R"], "*", ms=16, color=plotting.PALETTE[3], markeredgecolor=plotting.INK, label=f"optimum {ver['mass_kg']:.3f} kg, SF {ver['SF']:.2f}")
    ax.text(T_RANGE[0] + 0.15, R_RANGE[0] + 0.15, "white = infeasible", va="bottom", fontsize=8.5, color=plotting.INK2)
    ax.set_xlim(T_RANGE[0] - 0.2, T_RANGE[1] + 0.2)
    ax.set_ylim(R_RANGE[0] - 0.4, R_RANGE[1] + 0.4)
    ax.set_xlabel("plate thickness t (mm)")
    ax.set_ylabel("inner fillet radius R (mm)")
    ax.set_title(f"Design space: mass; SF = {SF_REQ:g} boundary (red), deflection = {DEFL_MAX:g} mm boundary (orange dashed)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, fontsize=8.5)
    plotting.save(fig, OUT / "design_space.png")

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    sc = ax.scatter(doe.mass_kg, doe.SF, c=doe.R, cmap=plotting.sequential_cmap(), s=42, edgecolor=plotting.INK2, linewidth=0.4)
    fig.colorbar(sc, ax=ax, label="fillet radius R (mm)")
    ax.axhline(SF_REQ, color=plotting.PALETTE[7], lw=1.5)
    ax.text(doe.mass_kg.min(), SF_REQ + 0.08, "required SF", va="bottom", fontsize=8.5, color=plotting.PALETTE[7])
    ax.plot(ver["mass_kg"], ver["SF"], "*", ms=15, color=plotting.PALETTE[3], markeredgecolor=plotting.INK)
    ax.plot(base["mass_kg"], base["SF"], "s", ms=9, color=plotting.INK)
    ax.set_xlabel("mass (kg)")
    ax.set_ylabel("safety factor on fillet stress")
    ax.set_title("Every FE run: a larger fillet buys safety factor without mass")
    plotting.save(fig, OUT / "mass_vs_sf.png")

    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    ax.plot(conv.nodes / 1e3, conv.sigma_fillet_MPa, marker="o", color=plotting.PALETTE[0], label="fillet (design metric)")
    ax.plot(conv.nodes / 1e3, conv.sigma_fixed_edge_MPa, marker="s", color=plotting.PALETTE[1], label="edge of fixed face (mesh-dependent)")
    ax.set_ylim(0, None)
    ax.set_xlabel("nodes (thousands)")
    ax.set_ylabel("max von Mises stress (MPa)")
    ax.set_title("Mesh refinement: fillet stress (design metric) vs stress at the fixed-face edge")
    ax.legend(loc="center left")
    plotting.save(fig, OUT / "stress_convergence.png")

    # ---------------- checks and summary
    c_fine, c_mid = conv.iloc[-1], conv.iloc[-2]
    checks = {
        "equilibrium: fixed-face reaction equals the applied load in every run (< 1e-6)": bool((doe.equilibrium_err < 1e-6).all() and conv.equilibrium_err.max() < 1e-6),
        "response surfaces fit the FE data (R2 > 0.99 for stress and deflection)": bool(r2_s > 0.99 and r2_d > 0.99),
        "verified optimum meets SF >= 2.0 and deflection <= limit in a direct FE run": bool(ver["SF"] >= SF_REQ and ver["deflection_mm"] <= DEFL_MAX),
        "verified optimum is lighter than the baseline design": bool(ver["mass_kg"] < base["mass_kg"]),
        "fillet stress converged: < 3 % change between the two finest meshes": bool(abs(c_fine.sigma_fillet_MPa / c_mid.sigma_fillet_MPa - 1) < 0.03),
    }
    fixed_edge_growth_pct = 100 * (c_fine.sigma_fixed_edge_MPa / conv.iloc[0].sigma_fixed_edge_MPa - 1)
    summary = {
        "requirements": dict(F_N=F, S_y_MPa=S_Y, SF_required=SF_REQ, deflection_max_mm=DEFL_MAX, t_range_mm=T_RANGE, R_range_mm=R_RANGE),
        "rsm": dict(stress_R2=r2_s, stress_max_resid=res_s, deflection_R2=r2_d, deflection_max_resid=res_d, mass_R2=r2_m),
        "surrogate_optimum": dict(t=t_opt, R=R_opt, mass_kg=float(best.fun)),
        "verified_optimum": ver, "trimmed": trimmed, "baseline": base,
        "mass_saving_pct": 100 * (1 - ver["mass_kg"] / base["mass_kg"]),
        "fixed_edge_stress_growth_pct_coarsest_to_finest": fixed_edge_growth_pct,
        "convergence": conv.to_dict(orient="records"), "checks": checks,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    doe[["t", "R", "elements", "mass_kg", "sigma_fillet_MPa", "SF", "deflection_mm", "equilibrium_err"]].to_markdown(OUT / "doe.md", index=False, floatfmt=".4g")
    conv[["size", "fillet_div", "nodes", "sigma_fillet_MPa", "sigma_fixed_edge_MPa", "deflection_mm"]].to_markdown(OUT / "convergence.md", index=False, floatfmt=".4g")
    print("\nChecks:")
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
