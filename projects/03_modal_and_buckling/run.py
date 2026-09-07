#!/usr/bin/env python
"""Project 03 - Modal and linear-buckling analysis of a steel cantilever, verified against closed forms.

Pipeline
    Gmsh structured hex mesh (C3D20R)  ->  CalculiX *FREQUENCY (10 modes, three meshes) and *BUCKLE
    (4 modes)  ->  every mode classified from its displacement pattern  ->  compared with the
    Euler-Bernoulli (and exact Timoshenko) cantilever frequencies, Saint-Venant torsion, the axial
    rod and Euler fixed-free column buckling  ->  mode-shape and convergence plots + PASS/FAIL checks.

Run:  python run.py        (one to two minutes on a laptop; results/ is rewritten, work/ is scratch)
"""
from __future__ import annotations

import json
import sys
import time
from math import cos, cosh, pi, sin, sinh, sqrt, tanh
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import NullFormatter
from scipy.optimize import brentq

from fealib import analytic as an, ccx, mesh as fm, plotting

HERE = Path(__file__).resolve().parent
WORK, OUT = HERE / "work", HERE / "results"

# ---------------------------------------------------------------- problem definition (SI: m, N, Pa, kg)
L, B, H = 0.400, 0.020, 0.010          # length, width (y), depth (z): the weak axis bends in z
E, NU, RHO = 200e9, 0.30, 7850.0       # structural steel
SEC = an.rect_section(B, H)            # A, I = B H^3 / 12 (weak axis), I_strong = H B^3 / 12
A, I_WEAK, I_STRONG = SEC["A"], SEC["I"], SEC["I_strong"]
G, KAPPA = an.shear_modulus(E, NU), an.shear_correction_rect(NU)
KGA = KAPPA * G * A                    # shear stiffness of the section (N), used by Timoshenko and Engesser
MESHES = [(40, 2, 1), (80, 4, 2), (160, 8, 4)]      # nx, ny, nz  ->  10, 5 and 2.5 mm cubic elements
N_MODES, N_BUCKLE = 10, 4
P_REF_PRESSURE = 1.0e7                 # reference pressure on the tip face (Pa), positive = compression
P_REF = P_REF_PRESSURE * B * H         # reference load (N): the expected first buckling factor is then ~2.6
SCALE_STUDY = [1.0, 1e2, 1e4, 1e6, 1e7, 1e8]        # reference pressures for the solver-conditioning study
SCALE_STUDY_MESH = MESHES[1]

WEAK, STRONG, TORSION, AXIAL, OTHER = "bending z (weak)", "bending y (strong)", "torsion x", "axial x", "other"
R_MAX = sqrt((B / 2) ** 2 + (H / 2) ** 2)          # corner distance from the section centroid


# ---------------------------------------------------------------- closed-form references
def torsion_constant_rect(a: float, b: float, terms: int = 50) -> float:
    """Saint-Venant torsion constant of an a x b rectangle (Timoshenko & Goodier series, a >= b)."""
    a, b = max(a, b), min(a, b)
    s = sum(tanh(n * pi * a / (2 * b)) / n ** 5 for n in range(1, 2 * terms, 2))
    return a * b ** 3 * (1.0 / 3.0 - 64.0 / pi ** 5 * (b / a) * s)


J_T = torsion_constant_rect(B, H)
I_P = A * (B ** 2 + H ** 2) / 12.0                 # polar second moment of area about the centroid


def torsion_frequencies(n: int) -> list[float]:
    """Fixed-free shaft, Saint-Venant torsion (free warping): f_n = (2n-1)/(4L) sqrt(G J / (rho I_p))."""
    return [(2 * k - 1) / (4.0 * L) * sqrt(G * J_T / (RHO * I_P)) for k in range(1, n + 1)]


def axial_frequencies(n: int) -> list[float]:
    """Fixed-free rod: f_n = (2n-1)/(4L) sqrt(E / rho)."""
    return [(2 * k - 1) / (4.0 * L) * sqrt(E / RHO) for k in range(1, n + 1)]


def _timoshenko_det(f: float, I: float) -> float:
    """Clamped-free Timoshenko beam frequency determinant (zero at a natural frequency).

    Non-dimensional form on xi = x/L: W'''' + (a + c) W'' - a (b - c) W = 0 with
    a = rho w^2 L^2 / (kappa G) (shear flexibility), b = kappa G A L^2 / (E I), c = rho I w^2 L^2 / (E I)
    (rotary inertia). W = C1 cosh(al xi) + C2 sinh(al xi) + C3 cos(be xi) + C4 sin(be xi); the rotation
    follows from psi' = W'' + a W. Clamped: W(0) = psi(0) = 0; free: psi'(1) = 0 and W'(1) - psi(1) = 0.
    """
    w = 2 * pi * f
    a = RHO * w ** 2 * L ** 2 / (KAPPA * G)
    b = KAPPA * G * A * L ** 2 / (E * I)
    c = RHO * I * w ** 2 * L ** 2 / (E * I)
    disc = sqrt((a + c) ** 2 + 4 * a * (b - c))
    al, be = sqrt((disc - (a + c)) / 2), sqrt((disc + (a + c)) / 2)
    kh, kt = (al ** 2 + a) / al, (be ** 2 - a) / be
    ch, sh, co, si = cosh(al), sinh(al), cos(be), sin(be)
    m11 = kh * al * ch + kt * be * co              # moment condition, coefficients of C1 and C2
    m12 = kh * (al * sh + be * si)
    m21 = -sh / al + si / be                       # shear condition (divided by a)
    m22 = -ch / al - (kh / kt) * co / be
    return (m11 * m22 - m12 * m21) / ch            # scaled to stay O(1)


def timoshenko_cantilever_frequencies(I: float, f_eb: list[float]) -> list[float]:
    """Exact Timoshenko frequencies: the root of the determinant just below each Euler-Bernoulli value."""
    out = []
    for fe in f_eb:
        grid = np.geomspace(0.80 * fe, 1.001 * fe, 600)
        vals = np.array([_timoshenko_det(f, I) for f in grid])
        idx = np.flatnonzero(np.sign(vals[:-1]) * np.sign(vals[1:]) < 0)
        out.append(brentq(_timoshenko_det, grid[idx[-1]], grid[idx[-1] + 1], args=(I,)))
    return out


F_EB_WEAK = an.cantilever_bending_frequencies(E, I_WEAK, RHO, A, L, 5)
F_EB_STRONG = an.cantilever_bending_frequencies(E, I_STRONG, RHO, A, L, 5)
F_TIMO = {WEAK: timoshenko_cantilever_frequencies(I_WEAK, F_EB_WEAK),
          STRONG: timoshenko_cantilever_frequencies(I_STRONG, F_EB_STRONG)}
F_THEORY = {WEAK: (F_EB_WEAK, "Euler-Bernoulli"), STRONG: (F_EB_STRONG, "Euler-Bernoulli"),
            TORSION: (torsion_frequencies(3), "Saint-Venant torsion"), AXIAL: (axial_frequencies(3), "axial rod")}
P_THEORY = {WEAK: an.fixed_free_buckling_loads(E, I_WEAK, L, N_BUCKLE),
            STRONG: [an.euler_critical_load(E, I_STRONG, L, K=2)] + an.fixed_free_buckling_loads(E, I_STRONG, L, N_BUCKLE)[1:]}


# ---------------------------------------------------------------- FE model
def build_mesh(nx: int, ny: int, nz: int) -> fm.Mesh:
    with fm.session("beam", order=2):
        return fm.box_hex(L, B, H, nx, ny, nz, order=2).with_element_type({"C3D20": "C3D20R"})


def deck_header(tag: str) -> str:
    return (f"** Project 03 - {tag}: steel cantilever {L * 1e3:g} x {B * 1e3:g} x {H * 1e3:g} mm, clamped on x = 0\n"
            "*INCLUDE, INPUT=mesh.inp\n"
            + ccx.material_block("STEEL", E=E, nu=NU, rho=RHO)
            + "*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL\n"
            + ccx.fmt_boundary("XMIN", 1, 3))


def solve_modal(m: fm.Mesh, tag: str):
    """*FREQUENCY step: returns (frequencies in Hz, one DISP block per mode, solver seconds)."""
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    m.write_inp(wd / "mesh.inp")
    (wd / "modal.inp").write_text(deck_header(tag) + f"*STEP\n*FREQUENCY\n{N_MODES}\n*NODE FILE\nU\n*END STEP\n")
    res = ccx.run_ccx(wd / "modal.inp")
    freqs = ccx.read_dat(res.dat).eigenfrequencies_hz
    blocks = ccx.read_frd(res.frd).all("DISP")
    assert len(freqs) == N_MODES == len(blocks), (len(freqs), len(blocks))
    assert all(abs(b.step_value - f) < 1e-3 * f for b, f in zip(blocks, freqs)), "frd/dat mode order mismatch"
    return freqs, blocks, res.seconds


def solve_buckle(m: fm.Mesh, tag: str, pressure: float):
    """*BUCKLE step with the compressive pressure applied inside the step: (factors, DISP blocks, seconds)."""
    wd = WORK / tag
    wd.mkdir(parents=True, exist_ok=True)
    m.write_inp(wd / "mesh.inp")
    (wd / "buckle.inp").write_text(deck_header(tag) + f"*STEP\n*BUCKLE\n{N_BUCKLE}\n"
                                   + m.dload_lines("XMAX", pressure) + "*NODE FILE\nU\n*END STEP\n")
    res = ccx.run_ccx(wd / "buckle.inp")
    factors = ccx.read_dat(res.dat).buckling_factors
    assert len(factors) == N_BUCKLE, factors
    blocks = ccx.read_frd(res.frd).all("DISP")[-N_BUCKLE:]        # the leading block is the static pre-load solution
    assert all(abs(b.step_value - lam) < 1e-3 * abs(lam) for b, lam in zip(blocks, factors)), "frd/dat mode order mismatch"
    return factors, blocks, res.seconds


# ---------------------------------------------------------------- mode post-processing
def station_profiles(m: fm.Mesh, block) -> pd.DataFrame:
    """Section-mean displacements and the section rotation about x at every x-station of the beam."""
    xyz = m.coords(block.node_ids)
    u = block.values[:, :3]
    ry, rz = xyz[:, 1] - B / 2, xyz[:, 2] - H / 2
    df = pd.DataFrame({"x": np.round(xyz[:, 0], 9), "ux": u[:, 0], "uy": u[:, 1], "uz": u[:, 2],
                       "mom": ry * u[:, 2] - rz * u[:, 1], "r2": ry ** 2 + rz ** 2})
    g = df.groupby("x", as_index=False).agg(ux=("ux", "mean"), uy=("uy", "mean"), uz=("uz", "mean"),
                                            mom=("mom", "sum"), r2=("r2", "sum"))
    g["theta"] = g.mom / g.r2                      # least-squares rigid rotation of the cross-section
    return g.drop(columns=["mom", "r2"])


def classify(prof: pd.DataFrame):
    """Dominant motion: section-mean ux / uy / uz amplitude, or the corner displacement due to twist."""
    amp = {AXIAL: prof.ux.abs().max(), STRONG: prof.uy.abs().max(), WEAK: prof.uz.abs().max(),
           TORSION: prof.theta.abs().max() * R_MAX}
    (top, a1), (_, a2) = sorted(amp.items(), key=lambda kv: kv[1], reverse=True)[:2]
    return (top if a1 > 2.0 * a2 else OTHER), amp


def analyse_modal(m: fm.Mesh, freqs, blocks):
    """Classify every mode and attach the closed-form reference for the k-th mode of its class."""
    rows, profiles, counts = [], [], {}
    for k, (f, blk) in enumerate(zip(freqs, blocks), start=1):
        prof = station_profiles(m, blk)
        label, _ = classify(prof)
        counts[label] = rank = counts.get(label, 0) + 1
        ref, name = F_THEORY.get(label, ([], ""))
        f_th = ref[rank - 1] if rank <= len(ref) else np.nan
        timo = F_TIMO.get(label, [])
        f_ti = timo[rank - 1] if rank <= len(timo) else np.nan
        rows.append(dict(mode=k, classification=label, rank=rank, f_FE_Hz=f, f_theory_Hz=f_th, theory=name,
                         err_pct=100 * (f - f_th) / f_th, f_timoshenko_Hz=f_ti,
                         err_vs_timoshenko_pct=100 * (f - f_ti) / f_ti))
        profiles.append(prof)
    return pd.DataFrame(rows), profiles


def analyse_buckling(m: fm.Mesh, factors, blocks, p_ref: float):
    rows, profiles, counts = [], [], {}
    for k, (lam, blk) in enumerate(zip(factors, blocks), start=1):
        prof = station_profiles(m, blk)
        label, _ = classify(prof)
        counts[label] = rank = counts.get(label, 0) + 1
        ref = P_THEORY.get(label, [])
        p_th = ref[rank - 1] if rank <= len(ref) else np.nan
        p_sh = p_th / (1.0 + p_th / KGA)               # Engesser: Euler load reduced by the shear flexibility
        p_fe = lam * p_ref
        rows.append(dict(mode=k, classification=label, rank=rank, buckling_factor=lam, P_cr_FE_N=p_fe,
                         P_cr_theory_N=p_th, err_pct=100 * (p_fe - p_th) / p_th,
                         P_cr_engesser_N=p_sh, err_vs_engesser_pct=100 * (p_fe - p_sh) / p_sh))
        profiles.append(prof)
    return pd.DataFrame(rows), profiles


def normalised_shape(prof: pd.DataFrame, comp: str, theory: np.ndarray) -> np.ndarray:
    """Deflection along the beam scaled to unit tip amplitude and sign-matched to the reference curve."""
    y = prof[comp].values / prof[comp].values[-1]
    return y if np.dot(y, theory) >= 0 else -y


def pick(table: pd.DataFrame, label: str, rank: int) -> pd.Series:
    hit = table[(table.classification == label) & (table["rank"] == rank)]
    if hit.empty:
        raise RuntimeError(f"no mode classified as {label!r} with rank {rank}")
    return hit.iloc[0]


def md(df: pd.DataFrame, path: Path) -> None:
    df = df.astype(object).where(df.notna(), None)          # NaN -> '-' in the markdown tables
    path.write_text(df.to_markdown(index=False, floatfmt=".6g", missingval="-") + "\n")


# ---------------------------------------------------------------- main
def main() -> int:
    t0 = time.perf_counter()
    plotting.setup()
    OUT.mkdir(exist_ok=True)
    print(f"theory: weak-axis EB {[round(f, 2) for f in F_EB_WEAK[:3]]} Hz, strong-axis EB {F_EB_STRONG[0]:.2f} Hz, "
          f"torsion {F_THEORY[TORSION][0][0]:.1f} Hz, axial {F_THEORY[AXIAL][0][0]:.1f} Hz; "
          f"P_cr weak {[round(p, 1) for p in P_THEORY[WEAK][:3]]} N, strong {P_THEORY[STRONG][0]:.1f} N")

    # ---------------- Part A + B on the mesh series
    runs = []
    for nx, ny, nz in MESHES:
        tag = f"hex_{nx}x{ny}x{nz}"
        m = build_mesh(nx, ny, nz)
        freqs, blocks, s_modal = solve_modal(m, tag)
        mtab, mprof = analyse_modal(m, freqs, blocks)
        factors, bblocks, s_buck = solve_buckle(m, tag, P_REF_PRESSURE)
        btab, bprof = analyse_buckling(m, factors, bblocks, P_REF)
        for t in (mtab, btab):
            t.insert(0, "mesh", f"{nx}x{ny}x{nz}")
            t.insert(1, "elements", m.n_elements)
            t.insert(2, "dofs", m.n_dofs())
        runs.append(dict(mesh=(nx, ny, nz), m=m, modal=mtab, modal_profiles=mprof, modal_s=s_modal,
                         buckle=btab, buckle_profiles=bprof, buckle_s=s_buck))
        print(f"{tag:14s} {m.n_elements:5d} C3D20R, {m.n_dofs():6d} dofs | modal {s_modal:5.1f}s: "
              + ", ".join(f"{r.f_FE_Hz:.1f} Hz {r.classification.split()[0]}-{r.classification.split()[1]}"
                          for r in mtab.itertuples())
              + f" | buckle {s_buck:5.1f}s: " + ", ".join(f"{r.P_cr_FE_N:.0f} N ({r.err_pct:+.2f}%)" for r in btab.itertuples()))
    fine = runs[-1]
    modal_all = pd.concat([r["modal"] for r in runs], ignore_index=True)
    buckle_all = pd.concat([r["buckle"] for r in runs], ignore_index=True)
    modal_all.to_csv(OUT / "modal_convergence.csv", index=False)
    buckle_all.to_csv(OUT / "buckling_convergence.csv", index=False)

    # ---------------- reference-load (solver conditioning) study for the buckling factor
    nx, ny, nz = SCALE_STUDY_MESH
    m_study = next(r["m"] for r in runs if r["mesh"] == SCALE_STUDY_MESH)
    study_rows = []
    for p in SCALE_STUDY:
        factors, blocks, s = solve_buckle(m_study, f"refload_{p:g}Pa_{nx}x{ny}x{nz}", p)
        label, _ = classify(station_profiles(m_study, blocks[0]))
        p_cr = factors[0] * p * B * H
        study_rows.append(dict(p_ref_Pa=p, P_ref_N=p * B * H, lambda_1=factors[0], mode1_classification=label,
                               P_cr1_FE_N=p_cr, err_vs_euler_pct=100 * (p_cr - P_THEORY[WEAK][0]) / P_THEORY[WEAK][0]))
        print(f"reference pressure {p:8.3g} Pa: lambda_1 = {factors[0]:.6g} ({label}) -> P_cr1 = {p_cr:.1f} N "
              f"({study_rows[-1]['err_vs_euler_pct']:+.2f}%)  {s:.1f}s")
    study = pd.DataFrame(study_rows)
    study.to_csv(OUT / "buckling_reference_load.csv", index=False)

    # ---------------- tables
    modal_fine = fine["modal"].drop(columns=["mesh", "elements", "dofs"]).copy()
    modal_fine["f_FE_Hz"] = modal_fine.f_FE_Hz.round(2)
    modal_fine["f_theory_Hz"] = modal_fine.f_theory_Hz.round(2)
    modal_fine["f_timoshenko_Hz"] = modal_fine.f_timoshenko_Hz.round(2)
    modal_fine["err_pct"] = modal_fine.err_pct.round(3)
    modal_fine["err_vs_timoshenko_pct"] = modal_fine.err_vs_timoshenko_pct.round(3)
    md(modal_fine.drop(columns=["rank"]), OUT / "modal_table.md")

    conv_rows = []
    for r in runs:
        t = r["modal"]
        row = dict(mesh=f"{r['mesh'][0]}x{r['mesh'][1]}x{r['mesh'][2]}", elements=r["m"].n_elements, dofs=r["m"].n_dofs())
        for label, rank, name in ((WEAK, 1, "weak1"), (WEAK, 2, "weak2"), (WEAK, 3, "weak3"), (STRONG, 1, "strong1")):
            s = pick(t, label, rank)
            row[f"f_{name}_Hz"] = round(s.f_FE_Hz, 3)
            row[f"err_{name}_pct"] = round(s.err_pct, 3)
        row["solve_s"] = round(r["modal_s"], 2)
        conv_rows.append(row)
    conv = pd.DataFrame(conv_rows)
    md(conv, OUT / "modal_convergence.md")

    buck_fine = fine["buckle"].drop(columns=["mesh", "elements", "dofs", "rank"]).copy()
    for c in ("P_cr_FE_N", "P_cr_theory_N", "P_cr_engesser_N"):
        buck_fine[c] = buck_fine[c].round(1)
    buck_fine["buckling_factor"] = buck_fine.buckling_factor.round(5)
    for c in ("err_pct", "err_vs_engesser_pct"):
        buck_fine[c] = buck_fine[c].round(3)
    md(buck_fine, OUT / "buckling_table.md")

    bconv_rows = []
    for r in runs:
        t = r["buckle"]
        row = dict(mesh=f"{r['mesh'][0]}x{r['mesh'][1]}x{r['mesh'][2]}", dofs=r["m"].n_dofs())
        for k in range(1, N_BUCKLE + 1):
            s = t[t["mode"] == k].iloc[0]
            row[f"P_cr{k}_N"] = round(s.P_cr_FE_N, 1)
            row[f"err{k}_pct"] = round(s.err_pct, 3)
        row["solve_s"] = round(r["buckle_s"], 2)
        bconv_rows.append(row)
    bconv = pd.DataFrame(bconv_rows)
    md(bconv, OUT / "buckling_convergence.md")

    study_md = study.copy()
    study_md["P_cr1_FE_N"] = study_md.P_cr1_FE_N.round(1)
    study_md["err_vs_euler_pct"] = study_md.err_vs_euler_pct.round(3)
    md(study_md, OUT / "buckling_reference_load.md")

    # ---------------- figure 1: weak-axis mode shapes (finest mesh) vs the Euler-Bernoulli shapes
    x_fine = np.linspace(0, L, 300)
    shapes = pd.DataFrame({"x_m": fine["modal_profiles"][0].x.values})
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.4), sharey=True)
    for k, ax in enumerate(axes, start=1):
        s = pick(fine["modal"], WEAK, k)
        prof = fine["modal_profiles"][int(s["mode"]) - 1]
        phi_nodes = an.cantilever_mode_shape(prof.x.values, L, k)
        fe = normalised_shape(prof, "uz", phi_nodes)
        shapes[f"fe_mode{k}"], shapes[f"theory_mode{k}"] = fe, phi_nodes
        ax.plot(x_fine * 1e3, an.cantilever_mode_shape(x_fine, L, k), color=plotting.INK2, lw=1.5, ls="--", label="Euler-Bernoulli")
        ax.plot(prof.x.values * 1e3, fe, color=plotting.PALETTE[k - 1], lw=1.6, marker="o", ms=3, markevery=8,
                label="CalculiX C3D20R")
        ax.axhline(0, color=plotting.GRID, lw=0.8)
        ax.set_title(f"weak-axis mode {k} (FE mode {int(s['mode'])})\n"
                     f"{s.f_FE_Hz:.2f} Hz, Euler-Bernoulli {s.f_theory_Hz:.2f} Hz", fontsize=10)
        ax.set_xlabel("x from clamp (mm)")
        ax.legend(loc="lower left" if k == 1 else "best")
    axes[0].set_ylabel("uz / uz(tip)")
    fig.suptitle(f"Weak-axis bending mode shapes, {fine['m'].n_elements} C3D20R elements", y=1.10)
    plotting.save(fig, OUT / "mode_shapes.png")
    shapes.to_csv(OUT / "mode_shapes.csv", index=False)

    # ---------------- figure 2: frequency error vs DOFs
    fig, ax = plt.subplots(figsize=(8.8, 4.4))
    ax.axhline(0, color=plotting.INK2, lw=1, ls="--", label="Euler-Bernoulli")
    series = ((WEAK, 1, "weak-axis mode 1", "-"), (WEAK, 2, "weak-axis mode 2", "-"),
              (WEAK, 3, "weak-axis mode 3", "-"), (STRONG, 1, "strong-axis mode 1", "--"))
    for i, (label, rank, name, ls) in enumerate(series):
        errs = [pick(r["modal"], label, rank).err_pct for r in runs]
        ax.plot(conv.dofs, errs, marker="o", ls=ls, color=plotting.PALETTE[i], label=name)
        f_ti, f_eb = F_TIMO[label][rank - 1], F_THEORY[label][0][rank - 1]
        ax.axhline(100 * (f_ti - f_eb) / f_eb, color=plotting.PALETTE[i], lw=1, ls=":", alpha=0.9)
    ax.plot([], [], color=plotting.INK2, lw=1, ls=":", label="Timoshenko - EB\n(shear + rotary inertia)")
    ax.set_xscale("log")
    ax.set_xticks(conv.dofs.tolist())
    ax.set_xticklabels([f"{d:,}\n({mm})" for d, mm in zip(conv.dofs, conv.mesh)])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("degrees of freedom (mesh nx x ny x nz)")
    ax.set_ylabel("frequency error vs Euler-Bernoulli (%)")
    ax.set_title("Modal convergence with mesh refinement (C3D20R)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    plotting.save(fig, OUT / "modal_convergence.png")

    # ---------------- figure 3: buckling mode shapes (finest mesh)
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.0), sharex=True)
    bshapes = pd.DataFrame({"x_m": fine["buckle_profiles"][0].x.values})
    for ax, row in zip(axes.ravel(), fine["buckle"].itertuples()):
        prof = fine["buckle_profiles"][row.mode - 1]
        comp = "uy" if row.classification == STRONG else "uz"
        n = row.rank
        th_nodes = 1 - np.cos((2 * n - 1) * pi * prof.x.values / (2 * L))
        fe = normalised_shape(prof, comp, th_nodes)
        bshapes[f"fe_mode{row.mode}_{comp}"], bshapes[f"theory_mode{row.mode}"] = fe, th_nodes
        ax.plot(x_fine * 1e3, 1 - np.cos((2 * n - 1) * pi * x_fine / (2 * L)), color=plotting.INK2, lw=1.5, ls="--",
                label=f"1 - cos({2 * n - 1} pi x / 2L)")
        ax.plot(prof.x.values * 1e3, fe, color=plotting.PALETTE[row.mode - 1], lw=1.6, marker="o", ms=3, markevery=8,
                label="CalculiX")
        ax.axhline(0, color=plotting.GRID, lw=0.8)
        ax.set_title(f"mode {row.mode}: {row.classification}\nP_cr = {row.P_cr_FE_N / 1e3:.2f} kN, "
                     f"theory {row.P_cr_theory_N / 1e3:.2f} kN ({row.err_pct:+.2f} %)", fontsize=10)
        ax.set_ylabel(f"{comp} / {comp}(tip)")
        ax.legend(loc="upper left", fontsize=8)
    for ax in axes[1]:
        ax.set_xlabel("x from clamp (mm)")
    fig.suptitle(f"Linear buckling mode shapes, fixed-free column, {fine['m'].n_elements} C3D20R elements", y=1.0)
    fig.tight_layout()
    plotting.save(fig, OUT / "buckling_modes.png")
    bshapes.to_csv(OUT / "buckling_modes.csv", index=False)

    # ---------------- objective checks
    weak3 = [pick(fine["modal"], WEAK, k) for k in (1, 2, 3)]
    strong1 = pick(fine["modal"], STRONG, 1)
    f1 = np.array([pick(r["modal"], WEAK, 1).f_FE_Hz for r in runs])
    d = np.diff(f1)
    b1 = fine["buckle"][fine["buckle"]["mode"] == 1].iloc[0]
    plateau = study[study.p_ref_Pa.between(1e4, 1e7)].P_cr1_FE_N
    checks = {
        "finest mesh: first three weak-axis bending frequencies within 1.5 % of Euler-Bernoulli":
            bool(all(abs(s.err_pct) < 1.5 for s in weak3)),
        "finest mesh: first strong-axis bending frequency within 1.5 % of Euler-Bernoulli":
            bool(abs(strong1.err_pct) < 1.5),
        "mode-1 frequency converges monotonically with mesh refinement (steps of one sign, shrinking)":
            bool((np.all(d < 0) or np.all(d > 0)) and abs(d[1]) < abs(d[0])),
        "mode-1 buckling load within 2 % of Euler fixed-free P_cr = pi^2 E I / (4 L^2)":
            bool(abs(b1.err_pct) < 2.0),
        "mode-1 buckling is about the weak axis (classified as z motion)":
            bool(b1.classification == WEAK),
        "buckling factor independent of the reference-load magnitude (10 kPa - 10 MPa agree within 0.1 %)":
            bool((plateau.max() - plateau.min()) / plateau.mean() < 1e-3),
        "all 10 modes on the finest mesh classified unambiguously (none 'other')":
            bool((fine["modal"].classification != OTHER).all()),
    }
    runtime = time.perf_counter() - t0
    summary = {
        "inputs": dict(L=L, B=B, H=H, E=E, nu=NU, rho=RHO, A=A, I_weak=I_WEAK, I_strong=I_STRONG, G=G, kappa=KAPPA,
                       J_torsion=J_T, I_polar=I_P, element="C3D20R", meshes=[list(mm) for mm in MESHES],
                       n_modes=N_MODES, n_buckling_modes=N_BUCKLE, reference_pressure_Pa=P_REF_PRESSURE, P_ref_N=P_REF),
        "theory": dict(f_weak_EB_Hz=F_EB_WEAK, f_strong_EB_Hz=F_EB_STRONG, f_weak_timoshenko_Hz=F_TIMO[WEAK],
                       f_strong_timoshenko_Hz=F_TIMO[STRONG], f_torsion_Hz=F_THEORY[TORSION][0],
                       f_axial_Hz=F_THEORY[AXIAL][0], P_cr_weak_N=P_THEORY[WEAK], P_cr_strong_N=P_THEORY[STRONG]),
        "finest_mesh": dict(mesh=list(fine["mesh"]), elements=fine["m"].n_elements, nodes=fine["m"].n_nodes,
                            dofs=fine["m"].n_dofs(), modal_solve_s=fine["modal_s"], buckle_solve_s=fine["buckle_s"]),
        "engineering": dict(axial_stress_at_Pcr1_MPa=b1.P_cr_FE_N / A / 1e6,
                            slenderness_ratio_weak=2 * L / sqrt(I_WEAK / A), kappa_G_A_N=KGA),
        "modal_finest": fine["modal"].to_dict("records"),
        "modal_convergence": conv.to_dict("records"),
        "buckling_finest": fine["buckle"].to_dict("records"),
        "buckling_convergence": bconv.to_dict("records"),
        "buckling_reference_load_study": study.to_dict("records"),
        "checks": checks,
        "runtime_s": runtime,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    print(f"\nTotal runtime {runtime:.1f} s\nChecks:")
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
