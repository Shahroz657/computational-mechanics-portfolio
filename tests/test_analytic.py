"""Closed-form helpers are the reference for every FE model, so they get their own tests."""
from math import pi, sqrt

import numpy as np
import pytest

from fealib import analytic as an


def test_cantilever_deflections_match_hand_values():
    sec = an.rect_section(0.02, 0.02)
    assert sec["I"] == pytest.approx(1.3333333e-8, rel=1e-6)
    assert an.cantilever_tip_deflection_eb(500, 0.4, 200e9, sec["I"]) == pytest.approx(4.0e-3, rel=1e-9)
    kappa = an.shear_correction_rect(0.3)
    assert kappa == pytest.approx(0.8497, abs=1e-4)
    d = an.cantilever_tip_deflection_timoshenko(500, 0.4, 200e9, sec["I"], an.shear_modulus(200e9, 0.3), sec["A"], kappa)
    assert d == pytest.approx(4.0076e-3, rel=1e-4)


def test_bending_stress_linear_in_x():
    sec = an.rect_section(0.02, 0.02)
    s = an.cantilever_bending_stress(np.array([0.0, 0.2, 0.4]), 500, 0.4, sec["c"], sec["I"])
    assert s[0] == pytest.approx(150e6, rel=1e-9)
    assert s[1] == pytest.approx(75e6, rel=1e-9)
    assert s[2] == pytest.approx(0.0, abs=1e-6)


def test_cantilever_modes_and_shapes():
    sec = an.rect_section(0.02, 0.01)
    f = an.cantilever_bending_frequencies(200e9, sec["I"], 7850, sec["A"], 0.4, n=3)
    f1_hand = (1.8751 ** 2 / (2 * pi * 0.4 ** 2)) * sqrt(200e9 * sec["I"] / (7850 * sec["A"]))
    assert f[0] == pytest.approx(f1_hand, rel=1e-4)
    assert f[1] / f[0] == pytest.approx((4.6941 / 1.8751) ** 2, rel=1e-4)
    x = np.linspace(0, 0.4, 5)
    for n in (1, 2, 3):
        phi = an.cantilever_mode_shape(x, 0.4, n)
        assert phi[0] == pytest.approx(0.0, abs=1e-12)
        assert phi[-1] == pytest.approx(1.0, rel=1e-12)


def test_euler_buckling_factors():
    assert an.euler_critical_load(200e9, 1e-9, 0.4, K=2.0) == pytest.approx(pi ** 2 * 200e9 * 1e-9 / 0.64)
    loads = an.fixed_free_buckling_loads(200e9, 1e-9, 0.4, n=3)
    assert loads[0] == pytest.approx(an.euler_critical_load(200e9, 1e-9, 0.4, K=2.0))
    assert loads[1] / loads[0] == pytest.approx(9.0)


def test_peterson_limits():
    assert float(an.kt_plate_hole_net(0.0)) == pytest.approx(3.0)
    assert float(an.kt_plate_hole_gross(0.0)) == pytest.approx(3.0)
    assert float(an.kirsch_sigma_yy_along_x(1.0, 1.0, 100.0)) == pytest.approx(300.0)
    assert float(an.kirsch_sigma_yy_along_x(1e6, 1.0, 100.0)) == pytest.approx(100.0, rel=1e-9)


def test_pin_fin_limits_and_balance():
    fin = an.pin_fin(h=50, k=200, D=0.006, L=0.06, T_base=100, T_inf=25)
    assert fin["temperature"](0.0) == pytest.approx(100.0)
    assert fin["T_tip"] < 100.0 and fin["T_tip"] > 25.0
    assert 0 < fin["eta_f"] <= 1.0
    long_fin = an.pin_fin(h=50, k=200, D=0.006, L=5.0, T_base=100, T_inf=25)
    m_inf = sqrt(50 * long_fin["P"] * 200 * long_fin["Ac"]) * 75           # infinite-fin heat rate
    assert long_fin["q_f"] == pytest.approx(m_inf, rel=1e-6)


def test_lame_boundary_conditions():
    sr, st = an.lame_stresses(np.array([0.05, 0.1]), 0.05, 0.1, 10e6)
    assert sr[0] == pytest.approx(-10e6, rel=1e-9)
    assert sr[1] == pytest.approx(0.0, abs=1.0)
