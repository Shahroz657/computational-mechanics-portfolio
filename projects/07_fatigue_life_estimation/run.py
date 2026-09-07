#!/usr/bin/env python
"""Project 07 - Fatigue life of the project-01 cantilever under a variable-amplitude tip load.

Pipeline
    synthetic 60 s tip-load history at 200 Hz  ->  root bending stress (0.3 MPa per N)
    ->  own ASTM E1049 rainflow count, cross-checked cycle-by-bin against the `rainflow` package
    ->  Shigley S-N line with Marin factors  ->  Goodman / Gerber / Soderberg / SWT mean-stress corrections
    ->  Palmgren-Miner damage per 60 s block and life  ->  plots + PASS/FAIL checks

Run:  python run.py        (a few seconds; results/ is rewritten)
Units: SI - lengths in m (mm where stated), loads in N, stresses in MPa, time in s.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rainflow as rf_pkg

from fealib import analytic as an, plotting
import fatigue as fat

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"

# ---------------------------------------------------------------- component (project 01 cantilever)
L, B, H = 0.400, 0.020, 0.020                       # m
SEC = an.rect_section(B, H)
MPA_PER_N = L * SEC["c"] / SEC["I"] / 1e6           # root bending stress per newton of tip load = 0.3 MPa/N

# ---------------------------------------------------------------- material (assumed: AISI 1045 normalised)
S_UT, S_Y = 585.0, 310.0                            # MPa
FINISH, LOADING, RELIABILITY = "machined", "bending", 0.90

# ---------------------------------------------------------------- load history
FS, T_BLOCK = 200.0, 60.0                           # Hz, s
P_MEAN = 150.0                                      # N
P1, F1 = 200.0, 4.0                                 # N, Hz   fast sinusoid
P2, F2 = 120.0, 0.3                                 # N, Hz   slow sinusoid
NOISE_SIGMA, MA_WINDOW, SEED = 40.0, 5, 7           # N, samples, rng seed
BIN = 10.0                                          # MPa, range-histogram bin for the rainflow comparison
METHODS = ("goodman", "gerber", "soderberg", "swt")
LABEL = {"goodman": "Goodman", "gerber": "Gerber", "soderberg": "Soderberg", "swt": "SWT"}

# ASTM E1049-85 section 5.4.4 worked example: history and its published rainflow counts (range: cycles)
ASTM_EXAMPLE = [-2.0, 1.0, -3.0, 5.0, -1.0, 3.0, -4.0, 4.0, -2.0]
ASTM_EXAMPLE_COUNTS = {3.0: 0.5, 4.0: 1.5, 6.0: 0.5, 8.0: 1.0, 9.0: 0.5}


def build_sn_curve() -> fat.SNCurve:
    k_a = fat.marin_surface(S_UT, FINISH)
    d_e = fat.equivalent_diameter_rect(B * 1e3, H * 1e3)
    k_b = fat.marin_size_bending(d_e)
    sn = fat.basquin_sn(S_UT, k_a=k_a, k_b=k_b, k_c=fat.marin_load(LOADING), k_d=1.0,
                        k_e=fat.marin_reliability(RELIABILITY))
    sn.marin["d_e_mm"] = d_e
    return sn


def build_history() -> pd.DataFrame:
    """Deterministic 60 s tip-load history: mean + two sinusoids + moving-average-filtered Gaussian noise."""
    t = np.arange(int(T_BLOCK * FS)) / FS
    rng = np.random.default_rng(SEED)
    noise = np.convolve(rng.normal(0.0, NOISE_SIGMA, t.size), np.ones(MA_WINDOW) / MA_WINDOW, mode="same")
    load = P_MEAN + P1 * np.sin(2 * np.pi * F1 * t) + P2 * np.sin(2 * np.pi * F2 * t) + noise
    return pd.DataFrame({"time_s": t, "load_N": load, "stress_MPa": MPA_PER_N * load})


def constant_amplitude_life(sn: fat.SNCurve, amplitude: float, infinite_life: bool, periods: int = 40) -> dict:
    """Rainflow + Miner life of a zero-mean cosine of the given amplitude sampled at 200 Hz, 4 Hz, starting and
    ending on a peak (so that every counted range is exactly twice the amplitude)."""
    t = np.linspace(0.0, periods / F1, int(periods * FS / F1) + 1)
    s = amplitude * np.cos(2 * np.pi * F1 * t)
    cyc = fat.rainflow_cycles(s)
    res = fat.miner_damage(cyc, sn, infinite_life=infinite_life, mean_stress="goodman", S_ut=S_UT, S_y=S_Y)
    N_ref = sn.N(amplitude, infinite_life=infinite_life)
    return dict(amplitude_MPa=amplitude, infinite_life=infinite_life, cycles_counted=res.cycles_per_block,
                N_curve=N_ref, life_cycles=res.life_cycles, rel_diff=abs(res.life_cycles / N_ref - 1.0))


def spectrum_by_bin(amplitudes: np.ndarray, counts: np.ndarray, width: float) -> tuple[np.ndarray, np.ndarray]:
    """Count per amplitude bin and the count-weighted mean amplitude of each occupied bin."""
    k = np.floor(amplitudes / width).astype(int)
    xs, ys = [], []
    for kk in np.unique(k):
        m = k == kk
        xs.append(counts[m].sum())
        ys.append(np.average(amplitudes[m], weights=counts[m]))
    return np.array(xs), np.array(ys)


def main() -> int:
    t0 = time.perf_counter()
    plotting.setup()
    OUT.mkdir(exist_ok=True)

    # ---------------- S-N line
    sn = build_sn_curve()
    S_1e3 = sn.S(sn.N_low)
    print(f"S-N: k_a={sn.marin['k_a']:.4f} k_b={sn.marin['k_b']:.4f} (d_e={sn.marin['d_e_mm']:.2f} mm) "
          f"k_c={sn.marin['k_c']} k_d={sn.marin['k_d']} k_e={sn.marin['k_e']}  S_e'={sn.S_e_prime:.1f}  "
          f"S_e={sn.S_e:.2f} MPa  f={sn.f:.4f}  a={sn.a:.2f} MPa  b={sn.b:.5f}")

    # ---------------- load history
    hist = build_history()
    hist.to_csv(OUT / "load_history.csv", index=False, float_format="%.6f")
    stress = hist.stress_MPa.to_numpy()
    print(f"history: {len(hist)} samples, load {hist.load_N.min():.1f}..{hist.load_N.max():.1f} N, "
          f"stress {stress.min():.2f}..{stress.max():.2f} MPa")

    # ---------------- rainflow: own implementation vs the rainflow package
    own = fat.rainflow_cycles(stress)
    ref = [(r, m, c) for r, m, c, _, _ in rf_pkg.extract_cycles(stress)]
    ranges, means, counts = fat.cycles_to_arrays(own)
    ranges_ref, _, counts_ref = fat.cycles_to_arrays(ref)
    n_bins = int(np.ceil(max(ranges.max(), ranges_ref.max()) / BIN))
    edges, h_own = fat.range_histogram(own, BIN, n_bins)
    _, h_ref = fat.range_histogram(ref, BIN, n_bins)
    total_own, total_ref = counts.sum(), counts_ref.sum()
    comp = pd.DataFrame({"range bin (MPa)": [f"{lo:.0f}-{hi:.0f}" for lo, hi in zip(edges[:-1], edges[1:])],
                         "own count": h_own, f"rainflow {rf_pkg.__version__} count": h_ref,
                         "difference": h_own - h_ref})
    comp.loc[len(comp)] = ["total", total_own, total_ref, total_own - total_ref]
    comp.to_csv(OUT / "rainflow_comparison.csv", index=False)
    (OUT / "rainflow_comparison.md").write_text(comp.to_markdown(index=False, floatfmt=".1f"))
    pd.DataFrame({"range_MPa": ranges, "mean_MPa": means, "count": counts}).to_csv(
        OUT / "rainflow_cycles.csv", index=False, float_format="%.6f")
    n_tp = fat.turning_points(stress)[0].size
    print(f"rainflow: {n_tp} reversals -> {len(own)} counted ranges, {total_own:.1f} cycles (own) vs "
          f"{total_ref:.1f} (package); max range {ranges.max():.2f} MPa; bins identical: {np.array_equal(h_own, h_ref)}")

    # ASTM E1049 worked example
    ex = fat.rainflow_cycles(ASTM_EXAMPLE)
    ex_counts: dict[float, float] = {}
    for r, _, c in ex:
        ex_counts[r] = ex_counts.get(r, 0.0) + c
    astm_ok = ex_counts == ASTM_EXAMPLE_COUNTS
    print(f"ASTM E1049 example: {ex_counts}  matches standard: {astm_ok}")

    # ---------------- Miner damage per 60 s block, four mean-stress corrections, with and without the cut-off
    dmg = {(m, il): fat.miner_damage(own, sn, infinite_life=il, mean_stress=m, S_ut=S_UT, S_y=S_Y)
           for m in METHODS for il in (True, False)}
    screen = pd.DataFrame([{
        "method": LABEL[m], "max sigma_ar (MPa)": dmg[(m, True)].sigma_ar_max,
        "S_e / max sigma_ar": sn.S_e / dmg[(m, True)].sigma_ar_max,
        "cycles above S_e": dmg[(m, True)].cycles_above_S_e,
        "D per block": dmg[(m, True)].damage,
        "life": "infinite" if not np.isfinite(dmg[(m, True)].life_blocks) else f"{dmg[(m, True)].life_blocks:.4g} blocks",
    } for m in METHODS])
    life = pd.DataFrame([{
        "method": LABEL[m], "D per block": dmg[(m, False)].damage,
        "life (blocks)": dmg[(m, False)].life_blocks, "life (h)": dmg[(m, False)].life_blocks * T_BLOCK / 3600.0,
        "life (cycles)": dmg[(m, False)].life_cycles,
        "life / Goodman life": dmg[(m, False)].life_blocks / dmg[("goodman", False)].life_blocks,
    } for m in METHODS])
    (OUT / "endurance_screening.md").write_text(screen.to_markdown(index=False, floatfmt=".4g"))
    (OUT / "miner_life.md").write_text(life.to_markdown(index=False, floatfmt=".4g"))
    life.to_csv(OUT / "miner_life.csv", index=False)
    lives = np.array([dmg[(m, False)].life_blocks for m in METHODS])
    spread = lives.max() / lives.min()
    print("\n" + screen.to_string(index=False) + "\n\n" + life.to_string(index=False) + f"\n\nspread max/min life = {spread:.3f}")

    # ---------------- constant-amplitude sanity checks
    ca = [constant_amplitude_life(sn, 200.0, infinite_life=False), constant_amplitude_life(sn, 300.0, infinite_life=True)]
    ca_df = pd.DataFrame(ca)
    (OUT / "constant_amplitude_check.md").write_text(ca_df.to_markdown(index=False, floatfmt=".6g"))
    print("\n" + ca_df.to_string(index=False))

    # ---------------- figure 1: history excerpt with reversals
    n5 = int(5 * FS)
    tp_idx, tp_val = fat.turning_points(stress)
    is_peak = np.empty(tp_val.size, dtype=bool)
    is_peak[1:] = tp_val[1:] > tp_val[:-1]
    is_peak[0] = tp_val[0] > tp_val[1]
    sel = tp_idx < n5
    fig, ax = plt.subplots(figsize=(7.6, 4.0))
    ax.plot(hist.time_s[:n5], stress[:n5], color=plotting.PALETTE[0], lw=1.2, label="root bending stress (0.3 MPa per N of tip load)")
    ax.plot(hist.time_s.values[tp_idx[sel & is_peak]], tp_val[sel & is_peak], ls="none", marker="^", ms=4.5,
            color=plotting.PALETTE[1], label="peaks")
    ax.plot(hist.time_s.values[tp_idx[sel & ~is_peak]], tp_val[sel & ~is_peak], ls="none", marker="v", ms=4.5,
            color=plotting.PALETTE[2], label="valleys")
    ax.axhline(MPA_PER_N * P_MEAN, color=plotting.INK2, lw=1, ls="--")
    ax.text(5.02, MPA_PER_N * P_MEAN, f"mean {MPA_PER_N * P_MEAN:.0f} MPa", va="center", fontsize=8.5, color=plotting.INK2)
    ax.set_xlim(0, 5)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("root bending stress (MPa)")
    ax.set_title("Load history excerpt (first 5 s of the 60 s block) with the reversals used by rainflow")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3)
    plotting.save(fig, OUT / "history_excerpt.png")

    # ---------------- figure 2: rainflow matrix (range vs mean)
    m_lo, m_hi = BIN * np.floor(means.min() / BIN), BIN * np.ceil(means.max() / BIN)
    edges_m = np.arange(m_lo, m_hi + BIN / 2, BIN)
    Hm, _, _ = np.histogram2d(means, ranges, bins=[edges_m, edges], weights=counts)
    Z = np.ma.masked_where(Hm.T == 0, Hm.T)
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    pc = ax.pcolormesh(edges_m, edges, Z, cmap=plotting.sequential_cmap(), edgecolors=plotting.GRID, linewidth=0.4)
    cb = fig.colorbar(pc, ax=ax, pad=0.02)
    cb.set_label("cycles per 60 s block")
    ax.set_xlabel("cycle mean stress (MPa)")
    ax.set_ylabel("cycle stress range (MPa)")
    ax.set_title(f"Rainflow matrix, {BIN:.0f} MPa bins ({total_own:.1f} cycles, own ASTM E1049 count)")
    ax.grid(False)
    plotting.save(fig, OUT / "rainflow_matrix.png")

    # ---------------- figure 3: S-N line and the block's cycle spectrum
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    N_line = np.logspace(3, 6, 60)
    ax.plot(N_line, sn.S(N_line), color=plotting.PALETTE[0], label="Basquin line  S = a N^b  (1e3 - 1e6 cycles)")
    ax.plot([1e6, 1e11], [sn.S_e, sn.S_e], color=plotting.PALETTE[0], ls="--", label=f"endurance limit S_e = {sn.S_e:.0f} MPa (cut-off)")
    N_ext = np.logspace(6, 11, 40)
    ax.plot(N_ext, sn.S(N_ext), color=plotting.PALETTE[0], ls=":", label="Basquin line extended below S_e (finite-life table)")
    x_raw, y_raw = spectrum_by_bin(0.5 * ranges, counts, BIN / 2)
    x_g, y_g = spectrum_by_bin(dmg[("goodman", False)].sigma_ar, counts, BIN / 2)
    ax.plot(x_raw, y_raw, ls="none", marker="o", mfc="none", mew=1.4, color=plotting.PALETTE[2], label="block cycles: raw amplitude vs count per block")
    ax.plot(x_g, y_g, ls="none", marker="o", color=plotting.PALETTE[1], label="block cycles: Goodman-equivalent amplitude")
    ax.plot([1e3], [S_1e3], marker="s", color=plotting.PALETTE[0], ls="none")
    ax.annotate(f"f S_ut = {S_1e3:.0f} MPa", (1e3, S_1e3), xytext=(6, 4), textcoords="offset points", fontsize=8.5, color=plotting.INK2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.3, 1e11)
    ax.set_ylim(10, 800)
    ax.set_xlabel("cycles  (markers: cycles per 60 s block;  line: cycles to failure)")
    ax.set_ylabel("stress amplitude (MPa)")
    ax.set_title("Corrected S-N line (Shigley) and the rainflow spectrum of one block")
    ax.legend(loc="lower left", fontsize=8)
    plotting.save(fig, OUT / "sn_curve.png")

    # ---------------- objective checks
    checks = {
        "own rainflow reproduces the ASTM E1049-85 worked example (ranges 3/4/6/8/9 -> 0.5/1.5/0.5/1.0/0.5 cycles)": bool(astm_ok),
        "rainflow total cycle count equals the rainflow package (|difference| < 1e-9)": bool(abs(total_own - total_ref) < 1e-9),
        "rainflow 10 MPa range histogram identical to the rainflow package in every bin": bool(np.array_equal(h_own, h_ref)),
        "constant-amplitude 200 MPa cosine: rainflow + Miner life equals N(200 MPa) within 1e-9 (Basquin line, no cut-off)": bool(ca[0]["rel_diff"] < 1e-9),
        "constant-amplitude 300 MPa cosine: rainflow + Miner life equals N(300 MPa) within 1e-9 (with the S_e cut-off)": bool(ca[1]["rel_diff"] < 1e-9),
        "endurance-limit screening: every corrected amplitude is below S_e, so the cut-off variant gives D = 0 for all methods": bool(
            all(dmg[(m, True)].damage == 0.0 and dmg[(m, True)].sigma_ar_max < sn.S_e for m in METHODS)),
        "Miner damage per block (Basquin line extended) is positive and finite for all four methods": bool(
            all(0 < dmg[(m, False)].damage < np.inf for m in METHODS)),
        "Goodman is at least as conservative as Gerber (life_goodman <= life_gerber)": bool(
            dmg[("goodman", False)].life_blocks <= dmg[("gerber", False)].life_blocks),
        "all four mean-stress methods agree within a factor of 10 (max/min life)": bool(spread < 10.0),
    }
    summary = {
        "inputs": dict(L=L, B=B, H=H, I=SEC["I"], MPa_per_N=MPA_PER_N, S_ut_MPa=S_UT, S_y_MPa=S_Y, finish=FINISH,
                       loading=LOADING, reliability=RELIABILITY, fs_Hz=FS, block_s=T_BLOCK, P_mean_N=P_MEAN,
                       P1_N=P1, f1_Hz=F1, P2_N=P2, f2_Hz=F2, noise_sigma_N=NOISE_SIGMA, ma_window=MA_WINDOW, seed=SEED),
        "sn_curve": dict(marin=sn.marin, S_e_prime_MPa=sn.S_e_prime, S_e_MPa=sn.S_e, f=sn.f, S_1e3_MPa=S_1e3,
                         a_MPa=sn.a, b=sn.b),
        "history": dict(samples=int(len(hist)), load_min_N=float(hist.load_N.min()), load_max_N=float(hist.load_N.max()),
                        stress_min_MPa=float(stress.min()), stress_max_MPa=float(stress.max())),
        "rainflow": dict(reversals=int(n_tp), counted_ranges=len(own), cycles_own=float(total_own),
                         cycles_package=float(total_ref), package_version=rf_pkg.__version__,
                         max_range_MPa=float(ranges.max()), bin_MPa=BIN, bins_identical=bool(np.array_equal(h_own, h_ref)),
                         astm_example_counts=ex_counts),
        "damage": {f"{m}_{'cutoff' if il else 'extended'}": dict(
            damage_per_block=dmg[(m, il)].damage, life_blocks=dmg[(m, il)].life_blocks,
            life_hours=dmg[(m, il)].life_blocks * T_BLOCK / 3600.0, life_cycles=dmg[(m, il)].life_cycles,
            sigma_ar_max_MPa=dmg[(m, il)].sigma_ar_max, cycles_above_S_e=dmg[(m, il)].cycles_above_S_e)
            for m in METHODS for il in (True, False)},
        "spread_max_over_min_life": float(spread),
        "constant_amplitude": ca,
        "runtime_s": time.perf_counter() - t0,
        "checks": checks,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=lambda o: float(o) if np.isscalar(o) else str(o)))
    print("\nChecks:")
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    print(f"\nruntime {summary['runtime_s']:.2f} s")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
