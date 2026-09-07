"""fatigue - stress-life (S-N) fatigue estimation for variable-amplitude loading.

Contents (SI throughout; stresses in MPa, lengths for the size factor in mm)
    turning_points, rainflow_cycles   ASTM E1049-85 rainflow cycle counting (section 5.4.4)
    marin_* helpers, basquin_sn       Shigley high-cycle S-N line with the Marin endurance-limit factors
    equivalent_fully_reversed         Goodman / Gerber / Soderberg / Smith-Watson-Topper mean-stress corrections
    miner_damage                      Palmgren-Miner linear damage summation per load block

References
    [1] ASTM E1049-85 (Reapproved 2017), Standard Practices for Cycle Counting in Fatigue Analysis, 5.4.4.
    [2] R. G. Budynas, J. K. Nisbett, Shigley's Mechanical Engineering Design, 9th ed., McGraw-Hill, 2011,
        Ch. 6: Sec. 6-7 (S_e'), 6-8 (S-N line, Fig. 6-18), 6-9 (Marin factors), 6-12 (mean stress), 6-15 (Miner).
    [3] N. E. Dowling, Mechanical Behavior of Materials, 4th ed., Pearson, 2013, Sec. 9.7 (mean stress, SWT).
    [4] K. N. Smith, P. Watson, T. H. Topper, "A stress-strain function for the fatigue of metals",
        Journal of Materials 5 (1970) 767-778.
    [5] M. A. Miner, "Cumulative damage in fatigue", J. Appl. Mech. 12 (1945) A159-A164.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import log10, sqrt

import numpy as np

# =============================================================================== rainflow counting
def turning_points(signal) -> tuple[np.ndarray, np.ndarray]:
    """Reduce a signal to its reversals (peaks and valleys) in the sense of ASTM E1049 [1].

    The first and last samples are always kept (start and end of the history), runs of equal values are
    collapsed to their first sample, and an interior sample is kept only where the slope changes sign.
    Returns ``(indices into signal, values)``.
    """
    x = np.asarray(signal, dtype=float).ravel()
    if x.size == 0:
        return np.array([], dtype=int), np.array([], dtype=float)
    keep = np.ones(x.size, dtype=bool)
    keep[1:] = x[1:] != x[:-1]                       # collapse plateaus
    idx, y = np.flatnonzero(keep), x[keep]
    if y.size < 3:
        return idx, y
    d = np.diff(y)
    interior = np.flatnonzero(d[:-1] * d[1:] < 0) + 1  # slope sign change
    sel = np.concatenate(([0], interior, [y.size - 1]))
    return idx[sel], y[sel]


def rainflow_cycles(signal) -> list[tuple[float, float, float]]:
    """Rainflow cycle counting per ASTM E1049-85 section 5.4.4 [1].

    The history is first reduced to peaks and valleys. Reversals are pushed one at a time onto a stack and,
    whenever at least three are present, the two most recent ranges are compared: X = |p[-1] - p[-2]|
    (newest) and Y = |p[-2] - p[-3]|.
      * X < Y            read the next reversal;
      * X >= Y and Y contains the starting point S (the oldest reversal on the stack): count Y as one
                         half cycle and discard S (the next point becomes the new starting point);
      * X >= Y otherwise: Y is a closed hysteresis loop - count it as one full cycle and remove its peak
                         and valley from the stack.
    When the data are exhausted every range still on the stack is counted as a half cycle.

    Returns a list of ``(range, mean, count)``: range = |peak - valley| (twice the amplitude),
    mean = (peak + valley) / 2, count = 1.0 for a full cycle and 0.5 for a half cycle.
    """
    _, pts = turning_points(signal)
    stack: list[float] = []
    cycles: list[tuple[float, float, float]] = []
    for p in pts:
        stack.append(float(p))
        while len(stack) >= 3:
            X = abs(stack[-1] - stack[-2])
            Y = abs(stack[-2] - stack[-3])
            if X < Y:
                break
            if len(stack) == 3:                       # Y contains the starting point
                cycles.append((Y, 0.5 * (stack[0] + stack[1]), 0.5))
                del stack[0]
            else:                                     # closed loop
                cycles.append((Y, 0.5 * (stack[-3] + stack[-2]), 1.0))
                del stack[-3:-1]
    for lo, hi in zip(stack[:-1], stack[1:]):         # residue: half cycles
        cycles.append((abs(hi - lo), 0.5 * (lo + hi), 0.5))
    return cycles


def cycles_to_arrays(cycles) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(ranges, means, counts)`` arrays from a list of (range, mean, count[, ...]) tuples."""
    if len(cycles) == 0:
        return np.array([]), np.array([]), np.array([])
    arr = np.array([c[:3] for c in cycles], dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2]


def range_histogram(cycles, bin_width: float, n_bins: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Cycle counts per range bin ``[k w, (k+1) w)``; returns ``(bin edges, counts)``."""
    ranges, _, counts = cycles_to_arrays(cycles)
    if n_bins is None:
        n_bins = int(np.ceil(ranges.max() / bin_width)) if ranges.size else 1
    edges = bin_width * np.arange(n_bins + 1)
    hist, _ = np.histogram(ranges, bins=edges, weights=counts)
    return edges, hist


# =============================================================================== S-N line (Shigley)
# Marin surface factor k_a = a S_ut^b, S_ut in MPa            [2, Table 6-2]
SURFACE_FINISH = {"ground": (1.58, -0.085), "machined": (4.51, -0.265), "cold-drawn": (4.51, -0.265),
                  "hot-rolled": (57.7, -0.718), "as-forged": (272.0, -0.995)}
# Marin load factor k_c                                        [2, Eq. 6-26]
LOAD_FACTOR = {"bending": 1.0, "axial": 0.85, "torsion": 0.59}
# Marin reliability factor k_e                                 [2, Table 6-5]
RELIABILITY_FACTOR = {0.50: 1.000, 0.90: 0.897, 0.95: 0.868, 0.99: 0.814, 0.999: 0.753,
                      0.9999: 0.702, 0.99999: 0.659, 0.999999: 0.620}


def marin_surface(S_ut: float, finish: str = "machined") -> float:
    """k_a = a S_ut^b with the coefficients of Shigley Table 6-2 (S_ut in MPa)."""
    a, b = SURFACE_FINISH[finish]
    return a * S_ut ** b


def equivalent_diameter_rect(b_mm: float, h_mm: float) -> float:
    """Equivalent rotating-beam diameter of a b x h rectangle in non-rotating bending,
    d_e = 0.808 sqrt(h b) (mm), Shigley Sec. 6-9 [2]."""
    return 0.808 * sqrt(h_mm * b_mm)


def marin_size_bending(d_mm: float) -> float:
    """k_b for bending or torsion, d in mm: 1.24 d^-0.107 (2.79 <= d <= 51), 1.51 d^-0.157 (51 < d <= 254),
    Shigley Eq. 6-20 [2]. Raises outside the tabulated range."""
    if 2.79 <= d_mm <= 51.0:
        return 1.24 * d_mm ** -0.107
    if 51.0 < d_mm <= 254.0:
        return 1.51 * d_mm ** -0.157
    raise ValueError(f"size factor: d = {d_mm:.2f} mm outside the 2.79-254 mm range of Shigley Eq. 6-20")


def marin_load(loading: str = "bending") -> float:
    """k_c: 1 bending, 0.85 axial, 0.59 torsion (Shigley Eq. 6-26)."""
    return LOAD_FACTOR[loading]


def marin_reliability(reliability: float = 0.50) -> float:
    """k_e from Shigley Table 6-5 (50 % -> 1.000, 90 % -> 0.897, 95 % -> 0.868, 99 % -> 0.814, ...)."""
    for r, k in RELIABILITY_FACTOR.items():
        if abs(r - reliability) < 1e-9:
            return k
    raise ValueError(f"reliability {reliability} not tabulated; choose one of {list(RELIABILITY_FACTOR)}")


def unnotched_endurance_limit(S_ut: float) -> float:
    """Rotating-beam endurance limit estimate S_e' = 0.5 S_ut (S_ut <= 1400 MPa) or 700 MPa (Shigley Eq. 6-8)."""
    return 0.5 * S_ut if S_ut <= 1400.0 else 700.0


def fatigue_strength_fraction(S_ut: float) -> float:
    """Fatigue strength fraction f = S_f(1e3) / S_ut, Shigley Fig. 6-18 [2, Sec. 6-8].

    The figure is generated from Eqs. 6-10 to 6-12 with the uncorrected S_e' = 0.5 S_ut:
        sigma_F' = S_ut + 345 MPa,   b = -log10(sigma_F' / S_e') / log10(2 N_e),  N_e = 1e6,
        f = (sigma_F' / S_ut) (2 N_f)^b,  N_f = 1e3.
    For S_ut < 490 MPa (70 kpsi) the figure is not drawn and Shigley recommends f = 0.9, used here.
    """
    if S_ut < 490.0:
        return 0.9
    sigma_F = S_ut + 345.0
    S_e_prime = unnotched_endurance_limit(S_ut)
    b = -log10(sigma_F / S_e_prime) / log10(2.0 * 1e6)
    return sigma_F / S_ut * (2.0 * 1e3) ** b


@dataclass(frozen=True)
class SNCurve:
    """High-cycle S-N line S = a N^b (Basquin form) between N_low = 1e3 and N_e = 1e6 cycles [2, Sec. 6-8]."""
    S_ut: float
    S_e_prime: float          # uncorrected rotating-beam endurance limit
    S_e: float                # corrected endurance limit at N_e cycles
    f: float                  # fatigue strength fraction at N_low cycles
    a: float
    b: float
    marin: dict = field(default_factory=dict)
    N_low: float = 1e3
    N_e: float = 1e6

    def S(self, N):
        """Fully reversed stress amplitude surviving N cycles on the Basquin line: S = a N^b."""
        return self.a * np.asarray(N, dtype=float) ** self.b

    def N(self, S, infinite_life: bool = False):
        """Cycles to failure at fully reversed amplitude S: N = (S / a)^(1/b) (Shigley Eq. 6-16).

        With ``infinite_life=True`` amplitudes at or below S_e return inf (endurance-limit cut-off);
        by default the Basquin line is extrapolated below S_e. Non-positive amplitudes return inf.
        """
        scalar = np.ndim(S) == 0
        S = np.atleast_1d(np.asarray(S, dtype=float))
        N = np.full(S.shape, np.inf)
        pos = S > 0
        N[pos] = (S[pos] / self.a) ** (1.0 / self.b)
        if infinite_life:
            N[S <= self.S_e] = np.inf
        return float(N[0]) if scalar else N

    def __call__(self, S, infinite_life: bool = False):
        return self.N(S, infinite_life)


def basquin_sn(S_ut: float, S_e: float | None = None, *, f: float | None = None, k_a: float = 1.0,
               k_b: float = 1.0, k_c: float = 1.0, k_d: float = 1.0, k_e: float = 1.0,
               k_f: float = 1.0) -> SNCurve:
    """Shigley high-cycle S-N line for steel [2, Secs. 6-7 to 6-9]:

        S_e' = 0.5 S_ut  (S_ut <= 1400 MPa, else 700 MPa)                  Eq. 6-8
        S_e  = k_a k_b k_c k_d k_e k_f S_e'    (Marin equation)              Eq. 6-18
        f    = fatigue strength fraction at 1e3 cycles (Fig. 6-18; 0.9 below 490 MPa)
        a    = (f S_ut)^2 / S_e,   b = -log10(f S_ut / S_e) / 3              Eqs. 6-14, 6-15
        S(N) = a N^b   (1e3 <= N <= 1e6),   N(S) = (S / a)^(1/b)             Eqs. 6-13, 6-16

    k_a surface, k_b size, k_c load, k_d temperature, k_e reliability, k_f miscellaneous (each 1 by
    default). Pass ``S_e`` to use a corrected endurance limit directly instead of the Marin product, and
    ``f`` to override the Fig. 6-18 fraction. Returns an :class:`SNCurve` whose ``N(S)`` is the callable.
    """
    S_e_prime = unnotched_endurance_limit(S_ut)
    marin = dict(k_a=k_a, k_b=k_b, k_c=k_c, k_d=k_d, k_e=k_e, k_f=k_f)
    if S_e is None:
        S_e = S_e_prime * k_a * k_b * k_c * k_d * k_e * k_f
    if f is None:
        f = fatigue_strength_fraction(S_ut)
    a = (f * S_ut) ** 2 / S_e
    b = -log10(f * S_ut / S_e) / 3.0
    return SNCurve(S_ut=S_ut, S_e_prime=S_e_prime, S_e=S_e, f=f, a=a, b=b, marin=marin)


# =============================================================================== mean-stress corrections
MEAN_STRESS_METHODS = ("goodman", "gerber", "soderberg", "swt")


def equivalent_fully_reversed(sigma_a, sigma_m, S_ut: float, S_y: float | None = None, method: str = "goodman"):
    """Equivalent fully reversed (zero-mean) amplitude sigma_ar of a cycle with amplitude sigma_a and mean
    sigma_m, i.e. the amplitude that would cause the same damage at R = -1 [2, Sec. 6-12; 3, Sec. 9.7].

        goodman     sigma_a/S_e + sigma_m/S_ut = 1      ->  sigma_ar = sigma_a / (1 - sigma_m/S_ut)
        gerber      sigma_a/S_e + (sigma_m/S_ut)^2 = 1  ->  sigma_ar = sigma_a / (1 - (sigma_m/S_ut)^2)
        soderberg   sigma_a/S_e + sigma_m/S_y = 1       ->  sigma_ar = sigma_a / (1 - sigma_m/S_y)
        swt         sigma_ar = sqrt(sigma_max sigma_a),  sigma_max = sigma_m + sigma_a   (Smith-Watson-Topper [4])

    Compressive means: Goodman, Gerber and Soderberg are evaluated with sigma_m = max(sigma_m, 0) - a
    compressive mean earns no benefit (conservative; the failure lines are not extended into compression).
    SWT keeps its own definition, which credits a compressive mean through sigma_max and gives sigma_ar = 0
    (no damage) when sigma_max <= 0. A tensile mean at or beyond the strength in the denominator
    (sigma_m >= S_ut, or >= S_y for Soderberg) returns inf (static failure of that cycle).
    Inputs may be scalars or arrays (broadcast); the return type follows the inputs.
    """
    sa = np.asarray(sigma_a, dtype=float)
    sm = np.asarray(sigma_m, dtype=float)
    method = method.lower()
    with np.errstate(divide="ignore", invalid="ignore"):
        if method == "goodman":
            denom = 1.0 - np.maximum(sm, 0.0) / S_ut
            out = np.where(denom > 0, sa / denom, np.inf)
        elif method == "gerber":
            denom = 1.0 - (np.maximum(sm, 0.0) / S_ut) ** 2
            out = np.where(denom > 0, sa / denom, np.inf)
        elif method == "soderberg":
            if S_y is None:
                raise ValueError("soderberg needs the yield strength S_y")
            denom = 1.0 - np.maximum(sm, 0.0) / S_y
            out = np.where(denom > 0, sa / denom, np.inf)
        elif method == "swt":
            smax = sm + sa
            out = np.sqrt(np.maximum(smax, 0.0) * sa)
        else:
            raise ValueError(f"unknown mean-stress method {method!r}; choose one of {MEAN_STRESS_METHODS}")
    out = np.where(sa > 0, out, 0.0)             # zero-amplitude 'cycles' never damage
    return float(out) if out.ndim == 0 else out


# =============================================================================== Palmgren-Miner
@dataclass
class MinerResult:
    damage: float                 # D per block
    life_blocks: float            # 1 / D  (inf when D = 0)
    cycles_per_block: float       # sum of the rainflow counts
    cycles_above_S_e: float       # counts whose corrected amplitude exceeds S_e
    sigma_ar_max: float
    sigma_ar: np.ndarray = field(repr=False)
    N: np.ndarray = field(repr=False)
    d: np.ndarray = field(repr=False)

    @property
    def life_cycles(self) -> float:
        return self.life_blocks * self.cycles_per_block


def miner_damage(cycles, sn_N, *, S_e: float | None = None, infinite_life: bool = True,
                 mean_stress: str | None = None, S_ut: float | None = None,
                 S_y: float | None = None) -> MinerResult:
    """Palmgren-Miner damage of one block of counted cycles, D = sum_i n_i / N_i [2, Sec. 6-15; 5].

    cycles        list of (range, mean, count) as returned by :func:`rainflow_cycles`; amplitude = range / 2.
    sn_N          an :class:`SNCurve` or any callable N(S) giving cycles to failure at fully reversed
                  amplitude S (Basquin line, no cut-off).
    mean_stress   None (use the raw amplitude) or a method of :func:`equivalent_fully_reversed`; S_ut
                  (and S_y for Soderberg) are then required.
    S_e           endurance limit for the cut-off; taken from ``sn_N.S_e`` when omitted.
    infinite_life True (default): cycles whose corrected amplitude does not exceed S_e are taken as
                  non-damaging (N_i = inf), Shigley's infinite-life criterion. False: the Basquin line is
                  extrapolated below S_e (the 'elementary' Miner rule, conservative for spectra whose large
                  cycles gradually lower the endurance limit) so every cycle contributes.

    Life in blocks is 1 / D (inf when D = 0); ``life_cycles`` = life_blocks x cycles per block.
    """
    ranges, means, counts = cycles_to_arrays(cycles)
    sigma_a = 0.5 * ranges
    if mean_stress is None:
        sigma_ar = sigma_a.copy()
    else:
        if S_ut is None:
            raise ValueError("mean-stress correction needs S_ut")
        sigma_ar = np.asarray(equivalent_fully_reversed(sigma_a, means, S_ut, S_y, mean_stress), dtype=float)
    if S_e is None:
        S_e = getattr(sn_N, "S_e", None)
    if infinite_life and S_e is None:
        raise ValueError("infinite_life=True needs S_e (pass it, or pass an SNCurve)")
    N_fn = sn_N.N if isinstance(sn_N, SNCurve) else sn_N
    with np.errstate(divide="ignore", invalid="ignore"):
        N_i = np.asarray(N_fn(sigma_ar), dtype=float)
        if infinite_life:
            N_i = np.where(sigma_ar > S_e, N_i, np.inf)
        d = np.where(np.isfinite(N_i), counts / N_i, np.where(N_i == 0, np.inf, 0.0))
    D = float(d.sum()) if d.size else 0.0
    above = float(counts[sigma_ar > S_e].sum()) if (S_e is not None and counts.size) else float("nan")
    return MinerResult(damage=D, life_blocks=(1.0 / D if D > 0 else np.inf), cycles_per_block=float(counts.sum()),
                       cycles_above_S_e=above, sigma_ar_max=float(sigma_ar.max()) if sigma_ar.size else 0.0,
                       sigma_ar=sigma_ar, N=N_i, d=d)
