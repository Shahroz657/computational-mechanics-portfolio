"""End-to-end: mesh, deck, solve, parse; skipped when no ccx executable is available."""
import numpy as np
import pytest

from fealib import analytic as an, ccx, mesh as fm

try:
    CCX = ccx.find_ccx()
except FileNotFoundError:
    CCX = None


@pytest.mark.skipif(CCX is None, reason="CalculiX ccx not installed")
def test_cantilever_matches_timoshenko(tmp_path):
    L, b, h, E, nu, P = 0.4, 0.02, 0.02, 200e9, 0.3, 500.0
    with fm.session("beam", order=2):
        m = fm.box_hex(L, b, h, 16, 2, 2, order=2).with_element_type({"C3D20": "C3D20R"})
    m.write_inp(tmp_path / "mesh.inp")
    loads = m.traction_loads("XMAX", [0, 0, -P / (b * h)])
    deck = ("*INCLUDE, INPUT=mesh.inp\n" + ccx.material_block("S", E=E, nu=nu)
            + "*SOLID SECTION, ELSET=EALL, MATERIAL=S\n*STEP\n*STATIC\n" + ccx.fmt_boundary("XMIN", 1, 3)
            + ccx.fmt_cload(loads) + "*NODE FILE\nU, S\n*NODE PRINT, NSET=XMIN, TOTALS=ONLY\nRF\n*END STEP\n")
    (tmp_path / "beam.inp").write_text(deck)
    res = ccx.run_ccx(tmp_path / "beam.inp")
    frd, dat = ccx.read_frd(res.frd), ccx.read_dat(res.dat)
    sec = an.rect_section(b, h)
    u_ref = an.cantilever_tip_deflection_timoshenko(P, L, E, sec["I"], an.shear_modulus(E, nu), sec["A"], an.shear_correction_rect(nu))
    u_fe = -frd.get("DISP").at(m.nsets["XMAX"])[:, 2].mean()
    assert u_fe == pytest.approx(u_ref, rel=0.02)
    assert dat.totals[0]["values"][2] == pytest.approx(P, rel=1e-6)
