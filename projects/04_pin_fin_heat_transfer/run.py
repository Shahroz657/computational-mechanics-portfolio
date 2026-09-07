#!/usr/bin/env python
"""Project 04 - Pin fin: 3D steady conduction + convection in CalculiX, verified against 1D fin theory.

Pipeline
    Gmsh OCC cylinder, second-order tets (C3D10)  ->  CalculiX *HEAT TRANSFER, STEADY STATE with *FILM
    on the lateral surface and the tip  ->  parse .frd (NT, RFL) and .dat (base reaction flux)  ->
    compare with Incropera's convective-tip fin (Table 3.4, case A): mesh convergence, centre-line
    profile, length sweep, energy balance  ->  plots + PASS/FAIL checks.

Run:  python run.py        (a few minutes on a laptop; results/ is rewritten, work/ is scratch)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize

from fealib import analytic as an, ccx, mesh as fm, plotting

HERE = Path(__file__).resolve().parent
WORK, OUT = HERE / "work", HERE / "results"

# ---------------------------------------------------------------- problem definition (SI: m, W, K; temperatures in degC)
D, L0 = 0.006, 0.060                    # pin diameter, baseline length (axis along +x, base disc at x = 0)
K, RHO, CP = 200.0, 2700.0, 900.0       # aluminium; rho and cp play no role in a steady state but complete the material
T_B, T_INF, H = 100.0, 25.0, 50.0       # base temperature, ambient, film coefficient on the lateral surface and the tip
SIZES = {"D/2": D / 2, "D/4": D / 4, "D/8": D / 8}       # target element edge lengths, convergence study
CONVERGED = "D/8"
LENGTHS = [0.010, 0.020, 0.040, 0.060, 0.080, 0.100, 0.120, 0.140]     # length sweep (60 mm reuses the converged run)
BI = H * (D / 2) / K                    # Biot number of the cross-section, must be << 1 for 1D fin theory
ETYPE = "C3D10"                         # CalculiX takes the structural name in a *HEAT TRANSFER step (DC3D10 gives identical results)
GAIN_STEP, GAIN_LIMIT = 0.020, 0.05     # "diminishing returns": adding 20 mm of fin buys less than 5 % heat rate


def reference(L: float) -> dict:
    return an.pin_fin(H, K, D, L, T_B, T_INF)


def build_mesh(L: float, size: float) -> fm.Mesh:
    """OCC cylinder along x with physical surfaces BASE (x = 0), TIP (x = L), LATERAL and volume FIN; the embedded
    axis line gives nodes exactly on the centre line for the T(x) comparison."""
    with fm.session("fin", order=2):                          # session() meshes single-threaded: reproducible meshes
        return fm.cylinder_tet(L, D, size, order=2, face_names=("BASE", "TIP", "LATERAL"), volume_name="FIN", axis_line=True)


def solve(L: float, size: float, tag: str):
    """One steady heat-transfer run. Returns a results row, the centre-line profile and (mesh, T by node)."""
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    m = build_mesh(L, size)
    m.with_element_type({"C3D10": ETYPE}).write_inp(wd / "mesh.inp")
    deck = (
        f"** Pin fin {tag}: T = {T_B} C on BASE, film h = {H} W/m2K to T_inf = {T_INF} C on LATERAL and TIP\n"
        "*INCLUDE, INPUT=mesh.inp\n"
        + ccx.material_block("AL", k=K, rho=RHO, cp=CP)
        + "*SOLID SECTION, ELSET=EALL, MATERIAL=AL\n*STEP\n*HEAT TRANSFER, STEADY STATE\n"
        + ccx.fmt_boundary("BASE", 11, 11, T_B)
        + m.film_lines("LATERAL", T_INF, H) + m.film_lines("TIP", T_INF, H)
        + "*NODE FILE\nNT, RFL\n*NODE PRINT, NSET=BASE, TOTALS=ONLY\nRFL\n*END STEP\n"
    )
    (wd / "fin.inp").write_text(deck)
    res = ccx.run_ccx(wd / "fin.inp")
    frd, dat = ccx.read_frd(res.frd), ccx.read_dat(res.dat)
    nt = frd.get("NDTEMP")
    T = {int(n): float(v) for n, v in zip(nt.node_ids, nt.column("T"))}
    ref = reference(L)

    # fin heat rate = total reaction flux at the prescribed-temperature base nodes (.dat TOTALS=ONLY line reads
    # "total heat generation for set BASE"; CalculiX signs it positive when heat flows into the body)
    q_base_dat = next(t for t in dat.totals if t["set"] == "BASE")["values"][0]
    rb = frd.get("RFL")
    rfl = {int(n): float(v) for n, v in zip(rb.node_ids, rb.column("RFL"))}
    base = m.nsets["BASE"]
    q_base_frd = float(sum(rfl[n] for n in base))
    q_f = abs(q_base_dat)
    # convective loss recomputed from the nodal temperatures: f_i = integral of h (T - T_inf) N_i over the film
    # faces (consistent nodal film loads, quadratic curved faces, 7-point Gauss rule); q_conv = sum of f_i
    theta = {n: T[n] - T_INF for n in T}
    f = {}
    for s in ("LATERAL", "TIP"):
        for n, v in m.surface_nodal_integral(s, theta).items():
            f[n] = f.get(n, 0.0) + H * v
    q_conv = sum(f.values())
    area_fe = m.surface_integral("LATERAL") + m.surface_integral("TIP")
    # CalculiX's RFL is the conduction flux vector K_cond.T: at a free film node it equals -f_i, at the base-ring
    # nodes it is the Dirichlet reaction minus the film load of the adjacent lateral faces. The reaction proper
    # (heat supplied at the base) is therefore the printed total plus the base-ring film share.
    film_share_base = sum(f.get(n, 0.0) for n in base)
    q_reaction = q_base_dat + film_share_base
    base_set = set(base)
    rfl_residual = max(abs(rfl[n] + f[n]) for n in f if n not in base_set)

    tip = m.nsets["TIP"]
    t_tip = float(np.mean([T[n] for n in tip]))
    axis = m.select_nodes(lambda x, y, z: abs(y) < 1e-9 and abs(z) < 1e-9)
    prof = pd.DataFrame({"x": [m.nodes[n][0] for n in axis], "T_fe": [T[n] for n in axis]}).sort_values("x", ignore_index=True)
    prof["T_1d"] = ref["temperature"](prof.x.values)
    prof["dT"] = prof.T_fe - prof.T_1d
    ids = np.array(list(T))
    dev_all = np.array([T[n] for n in ids]) - ref["temperature"](np.array([m.nodes[n][0] for n in ids]))

    row = dict(tag=tag, L_mm=L * 1e3, size_mm=size * 1e3, elements=m.n_elements, nodes=m.n_nodes,
               q_f_W=q_f, q_f_analytic_W=ref["q_f"], q_f_err_pct=100 * (q_f / ref["q_f"] - 1),
               T_tip_C=t_tip, T_tip_analytic_C=ref["T_tip"], T_tip_err_C=t_tip - ref["T_tip"],
               eta_f=q_f / (H * ref["A_f"] * (T_B - T_INF)), eta_f_analytic=ref["eta_f"], mL=ref["mL"],
               q_conv_W=q_conv, energy_balance_err_pct=100 * (q_conv / q_f - 1),
               q_reaction_W=q_reaction, q_reaction_err_pct=100 * (q_reaction / ref["q_f"] - 1),
               base_ring_film_share_W=film_share_base, reaction_closure_pct=100 * (q_reaction / q_conv - 1),
               rfl_film_residual_W=rfl_residual, max_nodal_film_load_W=max(abs(v) for v in f.values()),
               rfl_total_dat_W=q_base_dat, rfl_total_frd_W=q_base_frd,
               area_fe_mm2=area_fe * 1e6, area_analytic_mm2=ref["A_f"] * 1e6,
               dT_axis_max_C=float(prof.dT.abs().max()), dT_all_max_C=float(np.abs(dev_all).max()),
               solve_s=res.seconds)
    return row, prof, (m, T)


def md_table(df: pd.DataFrame, spec) -> str:
    """Markdown table with readable headers and one number format per column: spec = [(column, header, fmt), ...]."""
    cols = [c for c, _, _ in spec]
    return df[cols].rename(columns={c: h for c, h, _ in spec}).to_markdown(index=False, floatfmt=[f for _, _, f in spec])


def diminishing_returns(q_of_L, lengths) -> float | None:
    """Smallest L in ``lengths`` for which q(L + GAIN_STEP) / q(L) - 1 < GAIN_LIMIT, or None."""
    for L in lengths:
        q0, q1 = q_of_L(L), q_of_L(L + GAIN_STEP)
        if q0 is not None and q1 is not None and q1 / q0 - 1 < GAIN_LIMIT:
            return float(L)
    return None


def main() -> int:
    t_start = time.time()
    plotting.setup()
    OUT.mkdir(exist_ok=True)
    ref0 = reference(L0)
    print(f"analytic baseline: m={ref0['m']:.4f} 1/m  mL={ref0['mL']:.4f}  q_f={ref0['q_f']:.4f} W  "
          f"eta_f={ref0['eta_f']:.4f}  T_tip={ref0['T_tip']:.3f} C  Bi={BI:.2e}")

    # ---------------- study 1: mesh convergence at the baseline length
    conv_rows, profiles, fields = [], {}, {}
    for label, size in SIZES.items():
        row, prof, fld = solve(L0, size, f"conv_{label.replace('/', '')}")
        row = {"size": label, **row}
        conv_rows.append(row)
        profiles[label], fields[label] = prof, fld
        print(f"size {label:4s} h={size*1e3:.3f} mm  elems={row['elements']:6d}  q_f={row['q_f_W']:.4f} W "
              f"({row['q_f_err_pct']:+.3f}%)  reaction={row['q_reaction_W']:.4f} W ({row['q_reaction_err_pct']:+.3f}%)  "
              f"T_tip={row['T_tip_C']:.3f} C ({row['T_tip_err_C']:+.4f})  balance={row['energy_balance_err_pct']:+.3f}%  "
              f"closure={row['reaction_closure_pct']:+.1e}%  max|dT|axis={row['dT_axis_max_C']:.4f}  {row['solve_s']:.1f}s")
    conv = pd.DataFrame(conv_rows)
    conv.to_csv(OUT / "convergence.csv", index=False)
    best = conv[conv["size"] == CONVERGED].iloc[0]
    profiles[CONVERGED].to_csv(OUT / "axial_profile.csv", index=False)

    # ---------------- study 3: length sweep at the converged element size
    sweep_rows = []
    for L in LENGTHS:
        if abs(L - L0) < 1e-12:
            row = {k: v for k, v in best.items() if k != "size"}          # identical inputs: reuse the converged run
        else:
            row, _, _ = solve(L, SIZES[CONVERGED], f"sweep_L{L*1e3:.0f}")
        sweep_rows.append(row)
        print(f"L={row['L_mm']:5.0f} mm  mL={row['mL']:.3f}  elems={row['elements']:6d}  q_f={row['q_f_W']:.4f} W "
              f"({row['q_f_err_pct']:+.3f}%)  eta_f={row['eta_f']:.4f} (analytic {row['eta_f_analytic']:.4f})  "
              f"balance={row['energy_balance_err_pct']:+.3f}%  {row['solve_s']:.1f}s")
    sweep = pd.DataFrame(sweep_rows).astype({"elements": int, "nodes": int})
    sweep.to_csv(OUT / "length_sweep.csv", index=False)

    # diminishing returns: from the analytic curve on a 0.5 mm grid, and from the FE sweep (20 mm pairs)
    grid = np.round(np.arange(0.010, 0.400, 0.0005), 6)
    L_star_analytic = diminishing_returns(lambda L: reference(L)["q_f"], grid)
    q_fe = {round(r.L_mm): r.q_f_W for r in sweep.itertuples()}
    L_star_fe = diminishing_returns(lambda L: q_fe.get(round(L * 1e3)), LENGTHS)
    gains = pd.DataFrame([dict(L_from_mm=L * 1e3, L_to_mm=(L + GAIN_STEP) * 1e3,
                               gain_fe_pct=100 * (q_fe[round((L + GAIN_STEP) * 1e3)] / q_fe[round(L * 1e3)] - 1),
                               gain_analytic_pct=100 * (reference(L + GAIN_STEP)["q_f"] / reference(L)["q_f"] - 1))
                          for L in LENGTHS if round((L + GAIN_STEP) * 1e3) in q_fe])
    print(f"diminishing returns (< {GAIN_LIMIT:.0%} for +{GAIN_STEP*1e3:.0f} mm): analytic L* = {L_star_analytic*1e3:.1f} mm, "
          f"FE sweep L* = {L_star_fe*1e3:.0f} mm" if L_star_fe else "FE sweep does not reach the diminishing-returns point")

    # ---------------- figure 1: mesh convergence of q_f and T_tip
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.1))
    ax = axes[0]
    ax.axhline(ref0["q_f"], color=plotting.INK2, lw=1.2, ls="--", label="1D fin theory")
    ax.plot(conv.elements, conv.q_f_W, marker="o", color=plotting.PALETTE[0], label="base RFL total as printed (K_cond T)")
    ax.plot(conv.elements, conv.q_reaction_W, marker="s", ls="--", color=plotting.PALETTE[2], label="base reaction (RFL + base-ring film share)")
    for r in conv.itertuples():
        ax.annotate(f"{r.size}\n{r.q_f_err_pct:+.2f} %", (r.elements, r.q_f_W), textcoords="offset points", xytext=(0, -26),
                    ha="center", fontsize=8, color=plotting.INK2)
    ax.set_xscale("log")
    ax.set_xlabel("elements (C3D10)")
    ax.set_ylabel("fin heat rate q_f (W)")
    ax.set_title("Base heat rate")
    ax.margins(x=0.2, y=0.5)
    ax.legend(loc="upper left", fontsize=8)
    ax = axes[1]
    ax.axhline(ref0["T_tip"], color=plotting.INK2, lw=1.2, ls="--", label="1D fin theory")
    ax.plot(conv.elements, conv.T_tip_C, marker="o", color=plotting.PALETTE[1], label="CalculiX, mean over the tip nodes")
    for r in conv.itertuples():
        ax.annotate(f"{r.size}\n{r.T_tip_err_C:+.3f} C", (r.elements, r.T_tip_C), textcoords="offset points", xytext=(0, -26),
                    ha="center", fontsize=8, color=plotting.INK2)
    ax.set_xscale("log")
    ax.set_xlabel("elements (C3D10)")
    ax.set_ylabel("mean tip temperature (degC)")
    ax.set_title("Tip temperature")
    ax.margins(x=0.2, y=0.5)
    ax.legend(loc="upper left", fontsize=8)
    fig.suptitle(f"Mesh convergence, L = {L0*1e3:.0f} mm pin fin: element size D/2 -> D/4 -> D/8", fontweight="semibold")
    plotting.save(fig, OUT / "convergence.png")

    # ---------------- figure 2: centre-line temperature profile on the converged mesh
    prof = profiles[CONVERGED]
    xs = np.linspace(0, L0, 300)
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.6), sharex=True, gridspec_kw={"height_ratios": [3, 1.4]})
    ax = axes[0]
    ax.plot(xs * 1e3, ref0["temperature"](xs), color=plotting.INK2, lw=1.5, ls="--", label="1D fin theory, convective tip")
    sub = prof.iloc[::3]
    ax.plot(sub.x * 1e3, sub.T_fe, "o", ms=5.5, mfc="none", mew=1.4, color=plotting.PALETTE[0],
            label=f"CalculiX C3D10, centre-line nodes (every 3rd of {len(prof)})")
    ax.axhline(T_INF, color=plotting.INK2, lw=0.8, ls=":")
    ax.text(L0 * 1e3, T_INF + 1, "ambient", ha="right", va="bottom", fontsize=8.5, color=plotting.INK2)
    ax.set_ylabel("temperature (degC)")
    ax.set_title(f"Axial temperature profile, L = {L0*1e3:.0f} mm, mesh {CONVERGED}")
    ax.legend(loc="upper right")
    ax = axes[1]
    ax.axhline(0, color=plotting.INK2, lw=0.8)
    ax.plot(prof.x * 1e3, prof.dT, "o-", ms=3.5, lw=1, color=plotting.PALETTE[1])
    ax.set_xlabel("x from base (mm)")
    ax.set_ylabel("T_FE - T_1D (degC)")
    ax.set_title(f"Deviation from 1D theory (max |dT| = {prof.dT.abs().max():.3f} degC)", fontsize=10)
    plotting.save(fig, OUT / "axial_profile.png")

    # ---------------- figure 3: length sweep, q_f and eta_f (separate panels, one y-axis each)
    Ls = np.linspace(0.005, LENGTHS[-1] * 1.08, 300)
    q_an = np.array([reference(L)["q_f"] for L in Ls])
    eta_an = np.array([reference(L)["eta_f"] for L in Ls])
    to_mL, from_mL = (lambda L_mm: L_mm * 1e-3 * ref0["m"]), (lambda mL: mL / ref0["m"] * 1e3)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
    ax = axes[0]
    ax.plot(Ls * 1e3, q_an, color=plotting.INK2, lw=1.5, ls="--", label="1D fin theory")
    ax.plot(sweep.L_mm, sweep.q_f_W, "o", color=plotting.PALETTE[0], label=f"CalculiX C3D10, mesh {CONVERGED}")
    if L_star_analytic:
        ax.axvline(L_star_analytic * 1e3, color=plotting.PALETTE[3], lw=1.2, ls=":")
        ax.text(L_star_analytic * 1e3 + 2, q_an.min() + 0.06 * (q_an.max() - q_an.min()),
                f"+{GAIN_STEP*1e3:.0f} mm adds < {GAIN_LIMIT:.0%}\nbeyond L = {L_star_analytic*1e3:.0f} mm",
                fontsize=8.5, color=plotting.INK2)
    ax.set_xlabel("fin length L (mm)")
    ax.set_ylabel("fin heat rate q_f (W)")
    ax.set_title("Heat rate vs length", pad=26)
    ax.legend(loc="upper left")
    ax.secondary_xaxis("top", functions=(to_mL, from_mL)).set_xlabel("mL")
    ax = axes[1]
    ax.plot(Ls * 1e3, eta_an, color=plotting.INK2, lw=1.5, ls="--", label="1D fin theory")
    ax.plot(sweep.L_mm, sweep.eta_f, "o", color=plotting.PALETTE[1], label=f"CalculiX C3D10, mesh {CONVERGED}")
    ax.set_xlabel("fin length L (mm)")
    ax.set_ylabel("fin efficiency eta_f")
    ax.set_title("Efficiency vs length", pad=26)
    ax.set_ylim(0, 1.05)
    ax.legend(loc="lower left")
    ax.secondary_xaxis("top", functions=(to_mL, from_mL)).set_xlabel("mL")
    fig.suptitle("Length sweep: FE vs analytic (D = 6 mm, k = 200 W/m-K, h = 50 W/m2-K)", fontweight="semibold", y=1.16)
    plotting.save(fig, OUT / "length_sweep.png")

    # ---------------- figure 4: 3D temperature field on the converged mesh (surface faces, orthographic view at true
    # proportions, painter's algorithm) + deviation of every node from the 1D profile
    m, T = fields[CONVERGED]
    cmap = LinearSegmentedColormap.from_list("seq_blue", plotting.SEQ_BLUE)
    tris, vals = [], []
    for s in ("LATERAL", "TIP", "BASE"):
        for et, row in m.boundary[s]:
            tris.append([m.nodes[n] for n in row[:3]])
            vals.append(np.mean([T[n] for n in row]))
    tris, vals = np.array(tris) * 1e3, np.array(vals)
    az, el = np.radians(-20.0), np.radians(25.0)            # turn the tip towards the viewer, then tilt to show the top
    Rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, np.cos(el), -np.sin(el)], [0, np.sin(el), np.cos(el)]])
    view = tris @ (Rx @ Rz).T                                # screen: horizontal = x', vertical = z', depth = y' (viewer at -y')
    order = np.argsort(-view[:, :, 1].mean(axis=1))          # far faces first
    norm = Normalize(vmin=float(vals.min()), vmax=T_B)
    fig, axes = plt.subplots(2, 1, figsize=(9.6, 6.8), gridspec_kw={"height_ratios": [1.0, 1.5]})
    ax = axes[0]
    ax.add_collection(PolyCollection(view[order][:, :, [0, 2]], facecolors=cmap(norm(vals[order])), edgecolors="none"))
    ax.set_aspect("equal")
    ax.set_xlim(-6, L0 * 1e3 * np.cos(az) + 12)
    ax.set_ylim(-17, 11)
    ax.set_axis_off()
    base_top = np.array([0, 0, D / 2 * 1e3]) @ (Rx @ Rz).T
    tip_bot = np.array([L0 * 1e3, 0, -D / 2 * 1e3]) @ (Rx @ Rz).T
    ax.annotate(f"base, {T_B:.0f} degC (prescribed)", (base_top[0], base_top[2]), xytext=(2, 8.5), fontsize=9, color=plotting.INK,
                arrowprops=dict(arrowstyle="-", color=plotting.INK2, lw=0.8))
    ax.annotate(f"tip, {best.T_tip_C:.1f} degC (1D theory {ref0['T_tip']:.1f})", (tip_bot[0], tip_bot[2]), xytext=(tip_bot[0] - 22, -15),
                fontsize=9, color=plotting.INK, arrowprops=dict(arrowstyle="-", color=plotting.INK2, lw=0.8))
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.9, pad=0.01, aspect=12)
    cb.set_label("temperature (degC)")
    ax.set_title(f"Surface temperature, converged mesh: {m.n_elements} C3D10 elements ({m.n_nodes} nodes), true proportions")
    ax = axes[1]
    ids = np.array(list(T))
    xs_n = np.array([m.nodes[n][0] for n in ids]) * 1e3
    r_n = np.array([np.hypot(m.nodes[n][1], m.nodes[n][2]) for n in ids]) / (D / 2)
    dev = np.array([T[n] for n in ids]) - ref0["temperature"](xs_n / 1e3)
    sc = ax.scatter(xs_n, dev, c=r_n, cmap=cmap, s=2, linewidths=0, rasterized=True)
    ax.axhline(0, color=plotting.INK2, lw=0.8)
    cb = fig.colorbar(sc, ax=ax, pad=0.01, aspect=12)
    cb.set_label("r / (D/2)")
    ax.set_xlabel("x from base (mm)")
    ax.set_ylabel("T_FE - T_1D(x) (degC)")
    ax.set_title(f"Every node vs the 1D profile: max |dT| = {best.dT_all_max_C:.3f} degC, axis-to-surface spread "
                 f"{dev.max() - dev.min():.3f} degC (Bi = {BI:.1e})", fontsize=10)
    plotting.save(fig, OUT / "temperature_3d.png")

    # ---------------- engineering-note numbers (from this run / the closed forms, stored so the README quotes nothing else)
    sigma = 5.670374419e-8
    T_s, T_a = T_B + 273.15, T_INF + 273.15
    notes = dict(Bi=BI, dT_section_fe_C=float(dev.max() - dev.min()),
                 dT_section_estimate_C=H * (T_B - T_INF) * (D / 2) / (2 * K),      # parabolic radial profile at the base temperature
                 h_rad_W_m2K={f"eps={eps}": eps * sigma * (T_s ** 2 + T_a ** 2) * (T_s + T_a) for eps in (0.05, 0.8)},
                 effectiveness_baseline=ref0["effectiveness"])
    print(f"notes: section spread FE {notes['dT_section_fe_C']:.4f} C (estimate {notes['dT_section_estimate_C']:.4f} C); "
          f"h_rad at the base temperature: " + ", ".join(f"{k} {v:.2f} W/m2K" for k, v in notes["h_rad_W_m2K"].items()))

    # ---------------- objective checks
    dq = np.diff(conv.q_f_W.values)                    # successive changes with refinement
    monotone = bool(np.all(np.sign(dq) == np.sign(dq[0])) and np.all(np.abs(dq[1:]) < np.abs(dq[:-1])))
    checks = {
        f"converged ({CONVERGED}) fin heat rate q_f within 1.5 % of the analytic convective-tip solution": bool(abs(best.q_f_err_pct) < 1.5),
        f"converged ({CONVERGED}) mean tip temperature within 0.5 C of analytic": bool(abs(best.T_tip_err_C) < 0.5),
        "centre-line profile T(x): max |T_FE - T_1D| < 0.5 C on the converged mesh": bool(best.dT_axis_max_C < 0.5),
        "q_f converges monotonically (one-signed, shrinking changes D/2 -> D/4 -> D/8)": monotone,
        "length sweep: FE fin efficiency decreases monotonically with L": bool(np.all(np.diff(sweep.eta_f.values) < 0)),
        "energy balance: printed base flux equals the recomputed convective loss within 2 % (converged mesh)": bool(abs(best.energy_balance_err_pct) < 2),
        "energy balance, exact form: base reaction (RFL + base-ring film share) equals the convective loss within 0.01 % (converged mesh)": bool(abs(best.reaction_closure_pct) < 0.01),
        "RFL interpretation: at every free film node RFL = -(nodal film load) within 1 % of the largest nodal load (converged mesh)": bool(
            best.rfl_film_residual_W < 0.01 * best.max_nodal_film_load_W),
        "base heat rate from the .dat TOTALS line matches the sum of nodal RFL in the .frd (all runs, < 1e-4)": bool(
            (np.abs(pd.concat([conv, sweep]).eval("rfl_total_frd_W / rfl_total_dat_W - 1")) < 1e-4).all()),
    }
    runtime = time.time() - t_start
    summary = {
        "inputs": dict(D=D, L_baseline=L0, k=K, rho=RHO, cp=CP, T_base=T_B, T_inf=T_INF, h=H, element_type=ETYPE,
                       sizes_mm={k: v * 1e3 for k, v in SIZES.items()}, converged_size=CONVERGED,
                       sweep_lengths_mm=[L * 1e3 for L in LENGTHS]),
        "reference_baseline": dict(Bi=BI, m_per_m=ref0["m"], mL=ref0["mL"], q_f_W=ref0["q_f"], eta_f=ref0["eta_f"],
                                   effectiveness=ref0["effectiveness"], T_tip_C=ref0["T_tip"], A_f_mm2=ref0["A_f"] * 1e6),
        "converged_run": {k: (v.item() if hasattr(v, "item") else v) for k, v in best.items()},
        "convergence": conv.to_dict(orient="records"),
        "length_sweep": sweep.to_dict(orient="records"),
        "diminishing_returns": dict(step_mm=GAIN_STEP * 1e3, limit_pct=GAIN_LIMIT * 100,
                                    L_star_analytic_mm=L_star_analytic * 1e3 if L_star_analytic else None,
                                    L_star_fe_sweep_mm=L_star_fe * 1e3 if L_star_fe else None,
                                    gains=gains.to_dict(orient="records")),
        "rfl_sign_convention": "CalculiX RFL TOTALS at the fixed-temperature base nodes came out positive, i.e. positive = heat "
                               "flowing from the base into the fin; q_f is reported as the magnitude",
        "rfl_meaning": "RFL is the conduction flux vector K_cond.T: zero at interior nodes, minus the consistent nodal film "
                       "load at free film nodes, and reaction minus film load at the base-ring nodes; q_reaction_W adds "
                       "the base-ring film share back (base_ring_film_share_W)",
        "engineering_notes": notes,
        "runtime_s": runtime,
        "checks": checks,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    (OUT / "convergence.md").write_text(md_table(conv, [
        ("size", "size", ""), ("size_mm", "h (mm)", ".2f"), ("elements", "elements", ".0f"), ("nodes", "nodes", ".0f"),
        ("q_f_W", "q_f RFL (W)", ".4f"), ("q_f_err_pct", "err (%)", "+.3f"),
        ("q_reaction_W", "q_f reaction (W)", ".4f"), ("q_reaction_err_pct", "err (%)", "+.3f"),
        ("T_tip_C", "T_tip (C)", ".3f"), ("T_tip_err_C", "dT_tip (C)", "+.4f"), ("dT_axis_max_C", "max abs dT axis (C)", ".4f"),
        ("energy_balance_err_pct", "q_conv / q_f RFL - 1 (%)", "+.3f"), ("solve_s", "ccx (s)", ".2f")]))
    (OUT / "energy_balance.md").write_text(md_table(conv, [
        ("size", "size", ""), ("rfl_total_dat_W", "RFL total, base (W)", ".4f"), ("base_ring_film_share_W", "base-ring film share (W)", ".4f"),
        ("q_reaction_W", "reaction (W)", ".4f"), ("q_conv_W", "convective loss (W)", ".4f"), ("reaction_closure_pct", "closure (%)", ".1e"),
        ("rfl_film_residual_W", "max abs(RFL_i + f_i), free nodes (W)", ".1e"), ("max_nodal_film_load_W", "max abs f_i (W)", ".1e"),
        ("area_fe_mm2", "A_f mesh (mm2)", ".2f"), ("area_analytic_mm2", "A_f exact (mm2)", ".2f")]))
    (OUT / "length_sweep.md").write_text(md_table(sweep, [
        ("L_mm", "L (mm)", ".0f"), ("mL", "mL", ".3f"), ("elements", "elements", ".0f"), ("q_f_W", "q_f FE (W)", ".4f"),
        ("q_f_analytic_W", "q_f analytic (W)", ".4f"), ("q_f_err_pct", "err (%)", "+.3f"), ("eta_f", "eta_f FE", ".4f"),
        ("eta_f_analytic", "eta_f analytic", ".4f"), ("T_tip_C", "T_tip FE (C)", ".2f"), ("T_tip_analytic_C", "T_tip analytic (C)", ".2f"),
        ("energy_balance_err_pct", "q_conv / q_f - 1 (%)", "+.3f")]))
    (OUT / "diminishing_returns.md").write_text(md_table(gains, [
        ("L_from_mm", "L from (mm)", ".0f"), ("L_to_mm", "L to (mm)", ".0f"), ("gain_fe_pct", "q_f gain FE (%)", ".2f"),
        ("gain_analytic_pct", "q_f gain analytic (%)", ".2f")]))
    print(f"\nruntime {runtime:.0f} s\nChecks:")
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
