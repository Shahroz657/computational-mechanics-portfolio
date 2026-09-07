#!/usr/bin/env python
"""Project 06 - Reviewing an AI-generated FEA solution.

The review is executable: this script (1) reruns the AI solution exactly as written, (2) refines its
mesh to expose which of its numbers are mesh-dependent, (3) builds the corrected reference model, and
(4) applies a fixed list of objective checks to the AI solution, writing PASS/FAIL with evidence.
Units in this project follow the AI script (mm, N) so that the unit error can be demonstrated.
Run:  python run.py    (about 2 minutes)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import ai_solution as ai
from fealib import analytic as an, ccx, mesh as fm, plotting

HERE = Path(__file__).resolve().parent
WORK, OUT = HERE / "work", HERE / "results"
L, B, H, F = ai.L, ai.B, ai.H, ai.F
E_MPA, NU = 200e3, 0.30                         # consistent mm-N-MPa units
S_Y, S_UT = 275.0, 430.0
SEC = an.rect_section(B, H)
SIGMA_ROOT = float(an.cantilever_bending_stress(0.0, F, L, SEC["c"], SEC["I"]))
SIGMA_2H = float(an.cantilever_bending_stress(2 * H, F, L, SEC["c"], SEC["I"]))
U_TIP = an.cantilever_tip_deflection_timoshenko(F, L, E_MPA, SEC["I"], an.shear_modulus(E_MPA, NU), SEC["A"], an.shear_correction_rect(NU))
SIZES = [10.0, 6.0, 4.0, 2.5]


def top_fibre_stress_at_2H(m: fm.Mesh, stress, halfwidth: float) -> float:
    ids = m.select_nodes(lambda x, y, z: abs(z - H) < 1e-6 and abs(x - 2 * H) <= halfwidth)
    return float(stress.at(ids)[:, 0].mean())


def refine_ai_model(as_submitted: dict) -> pd.DataFrame:
    """Refine the AI model; the coarsest row is the submitted run itself (3D meshing is not bit-reproducible)."""
    rows = []
    for size in SIZES:
        r = as_submitted if size == as_submitted["mesh_size_mm"] else ai.run(size, WORK / f"ai_s{size:g}")
        m, stress, vm = r["_objects"]["mesh"], r["_objects"]["stress"], r["_objects"]["vm"]
        vm_by = dict(zip(stress.node_ids.tolist(), vm.tolist()))
        clamp_edge = max(vm_by[n] for n in m.nsets["XMIN"])
        rows.append(dict(mesh_size_mm=size, elements=r["elements"], nodes=r["nodes"],
                         sigma_load_node_MPa=vm_by[r["load_node"]], sigma_clamp_face_max_MPa=clamp_edge,
                         sigma_reported_max_MPa=r["sigma_max_MPa"], sigma_2H_top_fibre_MPa=top_fibre_stress_at_2H(m, stress, size / 2),
                         deflection_mm=r["deflection_mm"], reaction_z_N=r["reaction_z_N"]))
        print(f"AI model {size:4.1f} mm: {r['elements']:6d} C3D4  max={r['sigma_max_MPa']:7.1f} MPa  load node={vm_by[r['load_node']]:7.1f}  "
              f"clamp face={clamp_edge:7.1f}  x=2H top fibre={rows[-1]['sigma_2H_top_fibre_MPa']:6.1f} (theory {SIGMA_2H:.1f})  defl={r['deflection_mm']:.5f} mm")
    return pd.DataFrame(rows)


def corrected_model(n: int = 4) -> dict:
    wd = WORK / f"corrected_n{n}"
    wd.mkdir(parents=True, exist_ok=True)
    with fm.session("beam", order=2):
        m = fm.box_hex(L, B, H, 20 * n, n, n, order=2).with_element_type({"C3D20": "C3D20R"})
    m.write_inp(wd / "mesh.inp")
    loads = m.traction_loads("XMAX", [0.0, 0.0, -F / (B * H)])
    deck = ("*INCLUDE, INPUT=mesh.inp\n" + ccx.material_block("STEEL", E=E_MPA, nu=NU)
            + "*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL\n*STEP\n*STATIC\n" + ccx.fmt_boundary("XMIN", 1, 3)
            + ccx.fmt_cload(loads) + "*NODE FILE\nU, S\n*NODE PRINT, NSET=XMIN, TOTALS=ONLY\nRF\n*END STEP\n")
    (wd / "beam.inp").write_text(deck)
    res = ccx.run_ccx(wd / "beam.inp")
    frd, dat = ccx.read_frd(res.frd), ccx.read_dat(res.dat)
    disp, stress = frd.get("DISP"), frd.get("STRESS")
    top = m.select_nodes(lambda x, y, z: abs(z - H) < 1e-6 and abs(x - 2 * H) < 1e-6)
    sxx_2h = float(stress.at(top)[:, 0].mean())
    rf = next(t for t in dat.totals if t["label"].startswith("force"))["values"]
    return dict(elements=m.n_elements, nodes=m.n_nodes, element_type="C3D20R", sigma_2H_top_fibre_MPa=sxx_2h,
                sigma_root_beam_theory_MPa=SIGMA_ROOT, deflection_mm=float(-disp.at(m.nsets["XMAX"])[:, 2].mean()),
                deflection_theory_mm=U_TIP, reaction_z_N=rf[2], safety_factor_yield=S_Y / SIGMA_ROOT)


def review_checks(ai_result: dict, refine: pd.DataFrame, ref: dict) -> list[dict]:
    """Objective checks applied to the AI solution. Each returns PASS/FAIL plus the evidence used."""
    coarse, fine = refine.iloc[0], refine.iloc[-1]
    ratio_E_sigma = ai.E / ai_result["sigma_max_MPa"]
    x, y, z = ai_result["location_max"]
    at_load = ai_result["node_max"] == ai_result["load_node"]
    at_clamp = abs(x) < 1e-6
    growth = fine.sigma_reported_max_MPa / coarse.sigma_reported_max_MPa - 1
    checks = [
        dict(check="C1 unit consistency: E / sigma_max between 1e2 and 1e5 for a metal at working stress",
             passed=bool(1e2 < ratio_E_sigma < 1e5), evidence=f"E = {ai.E:.3g} with stresses in MPa gives E/sigma = {ratio_E_sigma:.2e}; deflection {ai_result['deflection_mm']:.5f} mm vs {U_TIP:.3f} mm expected (factor {U_TIP/ai_result['deflection_mm']:.0f})"),
        dict(check="C2 global equilibrium: clamp reaction equals applied load within 1e-6",
             passed=bool(abs(ai_result["reaction_z_N"] - F) / F < 1e-6), evidence=f"reaction {ai_result['reaction_z_N']:.6g} N vs {F} N applied (necessary, not sufficient)"),
        dict(check="C3 mesh convergence of the reported quantity: < 5 % change under 4x refinement",
             passed=bool(abs(growth) < 0.05), evidence=f"reported max stress {coarse.sigma_reported_max_MPa:.1f} -> {fine.sigma_reported_max_MPa:.1f} MPa from {coarse.mesh_size_mm:g} to {fine.mesh_size_mm:g} mm elements ({100*growth:+.0f} %)"),
        dict(check="C4 reported maximum is not read at a point load or a constraint edge (singular locations)",
             passed=bool(not (at_load or at_clamp)), evidence=f"max at node {ai_result['node_max']} at x={x:.0f}, y={y:.0f}, z={z:.0f} mm; load node = {ai_result['load_node']} -> {'point-load node' if at_load else ('clamp face' if at_clamp else 'interior')}"),
        dict(check="C5 hand-calculation plausibility: sigma_max within 25 % of M c / I and deflection within 20 % of beam theory",
             passed=bool(abs(ai_result["sigma_max_MPa"] / SIGMA_ROOT - 1) < 0.25 and abs(ai_result["deflection_mm"] / U_TIP - 1) < 0.20),
             evidence=f"sigma_max {ai_result['sigma_max_MPa']:.1f} MPa vs M c / I = {SIGMA_ROOT:.1f} MPa; deflection {ai_result['deflection_mm']:.5f} mm vs {U_TIP:.3f} mm"),
        dict(check="C6 element formulation suited to bending: no linear tetrahedra (C3D4) in a bending-dominated part",
             passed=False, evidence=f"C3D4 used; on the same mesh sizes the top-fibre stress at x = 2H reaches only {fine.sigma_2H_top_fibre_MPa:.0f} MPa of the {SIGMA_2H:.0f} MPa expected (see refinement table)"),
        dict(check="C7 correct strength criterion: static safety factor of a ductile part is based on yield, not ultimate",
             passed=False, evidence=f"SF = S_ut / sigma = {ai_result['safety_factor']:.2f} reported; with S_y = {S_Y:.0f} MPa and the beam-theory root stress the factor is {S_Y/SIGMA_ROOT:.2f}"),
        dict(check="C8 load application: resultant applied as a distributed traction or coupling, not a single nodal force on a solid mesh",
             passed=False, evidence=f"single *CLOAD on node {ai_result['load_node']}; stress at that node grows {coarse.sigma_load_node_MPa:.0f} -> {fine.sigma_load_node_MPa:.0f} MPa with refinement (singular)"),
    ]
    return checks


def main() -> int:
    plotting.setup()
    OUT.mkdir(exist_ok=True)
    ai_result = ai.run(ai.MESH_SIZE, WORK / "ai_as_submitted")
    (OUT / "ai_report.md").write_text(ai.report(ai_result))
    refine = refine_ai_model(ai_result)
    refine.to_csv(OUT / "ai_model_refinement.csv", index=False)
    ref = corrected_model()
    print(f"corrected model: {ref['elements']} {ref['element_type']}  sigma(2H) = {ref['sigma_2H_top_fibre_MPa']:.1f} MPa (theory {SIGMA_2H:.1f}), "
          f"tip = {ref['deflection_mm']:.3f} mm (theory {U_TIP:.3f}), SF_yield = {ref['safety_factor_yield']:.2f}")
    checks = review_checks(ai_result, refine, ref)

    # ---------------- figure: what converges and what does not
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    h = refine.mesh_size_mm.values
    ax.plot(h, refine.sigma_load_node_MPa, marker="o", color=plotting.PALETTE[7], label="von Mises at the point-load node")
    ax.plot(h, refine.sigma_clamp_face_max_MPa, marker="s", color=plotting.PALETTE[1], label="max von Mises on the clamped face")
    ax.plot(h, refine.sigma_2H_top_fibre_MPa, marker="^", color=plotting.PALETTE[0], label="bending stress at x = 2H (C3D4)")
    ax.axhline(SIGMA_2H, color=plotting.INK2, lw=1.2, ls="--")
    ax.text(h.min(), SIGMA_2H + 3, f"M c / I at x = 2H = {SIGMA_2H:.0f} MPa (converged C3D20R: {ref['sigma_2H_top_fibre_MPa']:.0f} MPa)", ha="right", va="bottom", fontsize=8.5, color=plotting.INK2)
    ax.invert_xaxis()
    ax.set_xticks(h)
    ax.set_xticklabels([f"{v:g}" for v in h])
    ax.set_ylim(0, 160)
    ax.set_xlabel("element size (mm), refining to the right")
    ax.set_ylabel("stress (MPa)")
    ax.set_title("AI model under refinement: point-load and clamp-edge stresses keep rising")
    ax.legend(loc="lower right", fontsize=8.5)
    plotting.save(fig, OUT / "singularity.png")

    comp = pd.DataFrame([
        dict(quantity="max stress reported (MPa)", ai=f"{ai_result['sigma_max_MPa']:.1f} (at node {ai_result['node_max']})", corrected=f"{SIGMA_ROOT:.1f} nominal at the root (M c / I); FE at 2H {ref['sigma_2H_top_fibre_MPa']:.1f} vs theory {SIGMA_2H:.1f}"),
        dict(quantity="tip deflection (mm)", ai=f"{ai_result['deflection_mm']:.5f}", corrected=f"{ref['deflection_mm']:.3f} (Timoshenko {U_TIP:.3f})"),
        dict(quantity="safety factor", ai=f"{ai_result['safety_factor']:.2f} (S_ut / sigma)", corrected=f"{ref['safety_factor_yield']:.2f} (S_y / sigma_root)"),
        dict(quantity="elements", ai=f"{ai_result['elements']} C3D4, single mesh", corrected=f"{ref['elements']} C3D20R, convergence shown in project 01"),
        dict(quantity="load application", ai="single nodal force", corrected="uniform traction via consistent nodal loads"),
        dict(quantity="units", ai="mm geometry with E = 200e9 (Pa)", corrected="mm-N-MPa: E = 200 000 MPa"),
    ])
    comp.to_markdown(OUT / "comparison.md", index=False)
    refine.to_markdown(OUT / "ai_model_refinement.md", index=False, floatfmt=".4g")
    pd.DataFrame(checks).to_markdown(OUT / "review_checks.md", index=False)
    summary = dict(ai_as_submitted={k: v for k, v in ai_result.items() if not k.startswith("_")},
                   reference=ref, theory=dict(sigma_root_MPa=SIGMA_ROOT, sigma_2H_MPa=SIGMA_2H, u_tip_mm=U_TIP),
                   refinement=refine.to_dict(orient="records"), checks=checks,
                   n_failed=sum(not c["passed"] for c in checks))
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print("\nReview checks on the AI solution:")
    for c in checks:
        print(f"  [{'PASS' if c['passed'] else 'FAIL'}] {c['check']}\n         {c['evidence']}")
    # the script succeeds when the review itself is consistent: the reference matches theory
    ok = abs(ref["sigma_2H_top_fibre_MPa"] / SIGMA_2H - 1) < 0.02 and abs(ref["deflection_mm"] / U_TIP - 1) < 0.01
    print(f"\nreference model agrees with theory: {'yes' if ok else 'NO'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
