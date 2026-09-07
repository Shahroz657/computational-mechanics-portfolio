"""The .frd/.dat parsers are exercised on hand-written snippets in the exact CalculiX layout."""
import textwrap

import numpy as np
import pytest

from fealib import ccx

FRD = """    1C
    2C                                       2                                     1
 -1         1 0.00000E+00 0.00000E+00 0.00000E+00
 -1         2 1.00000E+00 0.00000E+00 0.00000E+00
 -3
    3C                                       1                                     1
 -1         1    1    0    1
 -2         1         2         2         2         2         2         2         2
 -3
  100CL  101 1.000000000           2                     0    1           1
 -4  DISP        4    1
 -5  D1          1    2    1    0
 -5  D2          1    2    2    0
 -5  D3          1    2    3    0
 -5  ALL         1    2    0    0    1ALL
 -1         1 0.00000E+00 0.00000E+00 0.00000E+00
 -1         2 1.00000E-03-2.00000E-03 3.00000E-03
 -3
  100CL  101 1.000000000           2                     0    1           1
 -4  STRESS      6    1
 -5  SXX         1    4    1    1
 -5  SYY         1    4    2    2
 -5  SZZ         1    4    3    3
 -5  SXY         1    4    1    2
 -5  SYZ         1    4    2    3
 -5  SZX         1    4    3    1
 -1         1 1.00000E+02 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00
 -1         2 1.00000E+02 1.00000E+02 1.00000E+02 0.00000E+00 0.00000E+00 0.00000E+00
 -3
 9999
"""

DAT = """
                        S T E P       1

 total force (fx,fy,fz) for set XMIN and time  0.1000000E+01

        1.378385E-08 -2.243318E-09  5.000000E+02

     E I G E N V A L U E   O U T P U T

 MODE NO    EIGENVALUE                       FREQUENCY
                                    REAL PART            IMAGINARY PART
                                (RAD/TIME)      (CYCLES/TIME)  (RAD/TIME)

      1   0.1054213E+07  0.1026748E+04  0.1634127E+03  0.0000000E+00
      2   0.4000000E+08  0.6324555E+04  0.1006584E+04  0.0000000E+00

     B U C K L I N G   F A C T O R   O U T P U T

 MODE NO       BUCKLING
               FACTOR

      1   0.1234568E+02
      2   0.9876543E+02
"""


def test_read_frd_blocks(tmp_path):
    p = tmp_path / "job.frd"
    p.write_text(FRD)
    res = ccx.read_frd(p)
    assert res.nodes[2] == (1.0, 0.0, 0.0)
    disp = res.get("DISP")
    assert disp.components == ["D1", "D2", "D3"]
    assert disp.at(2).tolist() == pytest.approx([1e-3, -2e-3, 3e-3])
    stress = res.get("STRESS")
    assert stress.column("SXX").tolist() == [100.0, 100.0]
    vm = ccx.von_mises(stress.values)
    assert vm[0] == pytest.approx(100.0)          # uniaxial
    assert vm[1] == pytest.approx(0.0, abs=1e-9)  # hydrostatic


def test_read_dat_sections(tmp_path):
    p = tmp_path / "job.dat"
    p.write_text(DAT)
    d = ccx.read_dat(p)
    assert d.totals[0]["set"] == "XMIN" and d.totals[0]["values"][2] == pytest.approx(500.0)
    assert d.eigenfrequencies_hz == pytest.approx([163.4127, 1006.584])
    assert d.buckling_factors == pytest.approx([12.34568, 98.76543])


def test_deck_formatting():
    txt = ccx.fmt_elements("C3D20", {7: list(range(1, 21))})
    lines = txt.splitlines()
    assert lines[0] == "*ELEMENT, TYPE=C3D20, ELSET=EALL"
    assert lines[1].endswith(",") and len(lines) == 4          # 21 ints wrapped at 10 per line
    assert ccx.fmt_boundary("XMIN", 1, 3) == "*BOUNDARY\nXMIN, 1, 3\n"
    assert "*CLOAD\n5, 3, -500\n" == ccx.fmt_cload({5: np.array([0.0, 0.0, -500.0])})
