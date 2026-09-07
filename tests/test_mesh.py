"""Gmsh -> CalculiX conversion: set naming, face matching and consistent nodal loads."""
import numpy as np
import pytest

from fealib import mesh as fm


@pytest.fixture(scope="module")
def hex_mesh():
    with fm.session("box", order=2):
        return fm.box_hex(0.4, 0.02, 0.03, 8, 2, 3, order=2)


def test_box_hex_counts_and_sets(hex_mesh):
    m = hex_mesh
    assert list(m.elements) == ["C3D20"] and m.n_elements == 8 * 2 * 3
    assert set(m.nsets) >= {"XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX", "BODY"}
    assert len(m.faces["XMAX"]) == 2 * 3 and len(m.faces["ZMAX"]) == 8 * 2
    assert m.surface_measure("XMAX") == pytest.approx(0.02 * 0.03)
    assert m.surface_measure("ZMAX") == pytest.approx(0.4 * 0.02)


def test_consistent_loads_have_exact_resultant(hex_mesh):
    loads = hex_mesh.traction_loads("XMAX", [0.0, 0.0, -1000.0])
    total = np.sum(list(loads.values()), axis=0)
    assert total == pytest.approx([0.0, 0.0, -1000.0 * 0.02 * 0.03])
    # serendipity faces: corner nodes carry negative load, midside nodes positive
    vals = np.array([v[2] for v in loads.values()])
    assert (vals > 0).any() and (vals < 0).any()


def test_c3d20_permutation_keeps_positive_volume(hex_mesh):
    """Corner ordering must give a right-handed hexahedron (positive Jacobian), otherwise ccx aborts."""
    m = hex_mesh
    for nn in list(m.elements["C3D20"].values())[:5]:
        p = np.array([m.nodes[i] for i in nn[:8]])
        v = np.dot(np.cross(p[1] - p[0], p[3] - p[0]), p[4] - p[0])
        assert v > 0


def test_tet_mesh_faces():
    with fm.session("tet", order=2):
        m = fm.box_tet(0.1, 0.02, 0.02, 0.01, order=2)
    assert "C3D10" in m.elements
    assert m.surface_measure("XMAX") == pytest.approx(0.02 * 0.02, rel=1e-6)
    assert all(1 <= k <= 4 for _, k in m.faces["XMAX"])


def test_cylinder_surface_integral_measures_curved_faces():
    """Gauss quadrature on curved quadratic faces recovers the true lateral area; corner-only measure does not."""
    from math import pi
    with fm.session("cyl", order=2):
        m = fm.cylinder_tet(0.06, 0.006, 0.0015, order=2)
    exact = pi * 0.006 * 0.06
    assert m.surface_integral("LATERAL") == pytest.approx(exact, rel=1e-3)
    assert m.surface_measure("LATERAL") < exact                      # flat corner triangles under-measure a cylinder
    assert m.surface_integral("XMIN") == pytest.approx(pi * 0.003 ** 2, rel=2e-3)   # flat end disc
