"""fealib.analytic - closed-form reference solutions used to verify the finite-element models.

Sources: Timoshenko & Gere (beams, buckling), Cowper 1966 (shear coefficient), Blevins (cantilever
modes), Pilkey & Peterson (stress concentration), Incropera & DeWitt (fins), Lame (thick cylinders).
"""
from __future__ import annotations

from math import cos, cosh, pi, sin, sinh, sqrt

import numpy as np

# ------------------------------------------------------------------ beams
def rect_section(b: float, h: float) -> dict:
    """Area and second moments of a b x h rectangle (h is the depth in the bending direction)."""
    return {"A": b * h, "I": b * h ** 3 / 12.0, "I_strong": h * b ** 3 / 12.0, "c": h / 2.0}


def shear_modulus(E: float, nu: float) -> float:
    return E / (2.0 * (1.0 + nu))


def shear_correction_rect(nu: float) -> float:
    """Cowper's shear coefficient for a rectangular section."""
    return 10.0 * (1.0 + nu) / (12.0 + 11.0 * nu)


def cantilever_tip_deflection_eb(P, L, E, I):
    """Euler-Bernoulli tip deflection of a cantilever with a tip load."""
    return P * L ** 3 / (3.0 * E * I)


def cantilever_tip_deflection_timoshenko(P, L, E, I, G, A, kappa):
    """Timoshenko tip deflection: bending term plus transverse shear term."""
    return P * L ** 3 / (3.0 * E * I) + P * L / (kappa * G * A)


def cantilever_deflection_curve(x, P, L, E, I):
    x = np.asarray(x, dtype=float)
    return P * x ** 2 * (3.0 * L - x) / (6.0 * E * I)


def cantilever_bending_stress(x, P, L, c, I):
    """Bending stress magnitude at fibre distance c, tip load P (M = P (L - x))."""
    x = np.asarray(x, dtype=float)
    return P * (L - x) * c / I


# ------------------------------------------------------------------ cantilever bending modes
BETA_L_CANTILEVER = (1.87510407, 4.69409113, 7.85475744, 10.99554073, 14.13716839)


def cantilever_bending_frequencies(E, I, rho, A, L, n: int = 5) -> list[float]:
    """Natural frequencies (Hz) of an Euler-Bernoulli cantilever, first n bending modes."""
    return [(bl ** 2 / (2.0 * pi * L ** 2)) * sqrt(E * I / (rho * A)) for bl in BETA_L_CANTILEVER[:n]]


def cantilever_mode_shape(x, L, n: int) -> np.ndarray:
    """Normalised (unit tip) mode shape of the n-th bending mode of a cantilever."""
    bl = BETA_L_CANTILEVER[n - 1]
    beta = bl / L
    sigma = (sinh(bl) - sin(bl)) / (cosh(bl) + cos(bl))
    x = np.asarray(x, dtype=float)
    phi = np.cosh(beta * x) - np.cos(beta * x) - sigma * (np.sinh(beta * x) - np.sin(beta * x))
    tip = cosh(bl) - cos(bl) - sigma * (sinh(bl) - sin(bl))
    return phi / tip


# ------------------------------------------------------------------ Euler buckling
def euler_critical_load(E, I, L, K: float = 1.0) -> float:
    """Euler critical load with effective-length factor K (pinned-pinned 1, fixed-free 2, fixed-fixed 0.5)."""
    return pi ** 2 * E * I / (K * L) ** 2


def fixed_free_buckling_loads(E, I, L, n: int = 3) -> list[float]:
    """Higher buckling loads of a fixed-free column: P_n = (2n-1)^2 pi^2 E I / (4 L^2)."""
    return [(2 * k - 1) ** 2 * pi ** 2 * E * I / (4.0 * L ** 2) for k in range(1, n + 1)]


# ------------------------------------------------------------------ stress concentration
def kt_plate_hole_net(d_over_W):
    """Peterson/Pilkey fit for a plate with a central circular hole in tension, referenced to the
    *net-section* nominal stress (valid for 0 <= d/W <= 1)."""
    r = np.asarray(d_over_W, dtype=float)
    return 3.000 - 3.140 * r + 3.667 * r ** 2 - 1.527 * r ** 3


def kt_plate_hole_gross(d_over_W):
    """Same factor referenced to the gross (far-field) stress: sigma_max / sigma_0."""
    r = np.asarray(d_over_W, dtype=float)
    return kt_plate_hole_net(r) / (1.0 - r)


def kirsch_sigma_yy_along_x(x, a, sigma0):
    """Kirsch: stress in the loading direction along the net section of an infinite plate, x >= a."""
    x = np.asarray(x, dtype=float)
    return sigma0 * (1.0 + a ** 2 / (2.0 * x ** 2) + 3.0 * a ** 4 / (2.0 * x ** 4))


# ------------------------------------------------------------------ fins (Incropera, convective tip)
def pin_fin(h, k, D, L, T_base, T_inf) -> dict:
    """Straight pin fin with convection from the tip (Incropera Table 3.4, case A)."""
    P = pi * D
    Ac = pi * D ** 2 / 4.0
    m = sqrt(h * P / (k * Ac))
    M = sqrt(h * P * k * Ac) * (T_base - T_inf)
    H = h / (m * k)
    denom = cosh(m * L) + H * sinh(m * L)
    q_f = M * (sinh(m * L) + H * cosh(m * L)) / denom
    A_f = P * L + Ac
    eta_f = q_f / (h * A_f * (T_base - T_inf))
    q_no_fin = h * Ac * (T_base - T_inf)

    def temperature(x):
        x = np.asarray(x, dtype=float)
        return T_inf + (T_base - T_inf) * (np.cosh(m * (L - x)) + H * np.sinh(m * (L - x))) / denom

    return {"m": m, "mL": m * L, "q_f": q_f, "eta_f": eta_f, "effectiveness": q_f / q_no_fin,
            "T_tip": float(temperature(L)), "temperature": temperature, "A_f": A_f, "Ac": Ac, "P": P}


# ------------------------------------------------------------------ thick-walled cylinder (Lame)
def lame_stresses(r, a, b, p_i, p_o=0.0):
    """Radial and hoop stress in a thick cylinder, inner radius a, outer b, pressures p_i, p_o."""
    r = np.asarray(r, dtype=float)
    A = (p_i * a ** 2 - p_o * b ** 2) / (b ** 2 - a ** 2)
    B = (p_i - p_o) * a ** 2 * b ** 2 / (b ** 2 - a ** 2)
    return A - B / r ** 2, A + B / r ** 2
