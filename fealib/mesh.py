"""fealib.mesh - Gmsh helpers that produce CalculiX-ready meshes.

A :class:`Mesh` carries nodes, elements grouped by CalculiX type, node/element sets named after the
Gmsh physical groups, and element-face lists for every physical *surface* (needed for pressure,
film and consistent traction loads). Node ordering is permuted from Gmsh to CalculiX convention.
"""
from __future__ import annotations

import contextlib
import os
from dataclasses import dataclass, field, replace
from math import sqrt
from pathlib import Path

import gmsh
import numpy as np

from . import ccx as _ccx

# gmsh element type id -> (CalculiX type, permutation gmsh -> ccx node order)
GMSH_TO_CCX = {
    2: ("CPS3", [0, 1, 2]),
    3: ("CPS4", [0, 1, 2, 3]),
    9: ("CPS6", [0, 1, 2, 3, 4, 5]),
    16: ("CPS8", [0, 1, 2, 3, 4, 5, 6, 7]),
    4: ("C3D4", [0, 1, 2, 3]),
    11: ("C3D10", [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]),
    5: ("C3D8", [0, 1, 2, 3, 4, 5, 6, 7]),
    17: ("C3D20", [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15]),
    6: ("C3D6", [0, 1, 2, 3, 4, 5]),
    18: ("C3D15", [0, 1, 2, 3, 4, 5, 6, 9, 7, 12, 14, 13, 8, 10, 11]),
}

# CalculiX face numbering (0-based corner indices, in F1, F2, ... order)
_HEX = [(0, 1, 2, 3), (4, 7, 6, 5), (0, 4, 5, 1), (1, 5, 6, 2), (2, 6, 7, 3), (3, 7, 4, 0)]
_TET = [(0, 1, 2), (0, 3, 1), (1, 3, 2), (2, 3, 0)]
_WEDGE = [(0, 1, 2), (3, 5, 4), (0, 3, 4, 1), (1, 4, 5, 2), (2, 5, 3, 0)]
_TRI = [(0, 1), (1, 2), (2, 0)]
_QUAD = [(0, 1), (1, 2), (2, 3), (3, 0)]
CCX_FACES = {
    "C3D8": _HEX, "C3D8R": _HEX, "C3D8I": _HEX, "C3D20": _HEX, "C3D20R": _HEX,
    "C3D4": _TET, "C3D10": _TET, "C3D6": _WEDGE, "C3D15": _WEDGE,
    "CPS3": _TRI, "CPS6": _TRI, "CPE3": _TRI, "CPE6": _TRI, "CAX3": _TRI, "CAX6": _TRI,
    "CPS4": _QUAD, "CPS8": _QUAD, "CPS4R": _QUAD, "CPS8R": _QUAD,
    "CPE4": _QUAD, "CPE8": _QUAD, "CAX4": _QUAD, "CAX8": _QUAD,
}
# thermal element names share the structural face numbering
for _k in list(CCX_FACES):
    if _k.startswith("C3D"):
        CCX_FACES["D" + _k] = CCX_FACES[_k]

# consistent nodal fractions of a uniform traction resultant, per gmsh boundary element type
_LUMP = {
    1: [0.5, 0.5],                                   # 2-node line
    8: [1 / 6, 1 / 6, 2 / 3],                        # 3-node line
    2: [1 / 3] * 3,                                  # 3-node triangle
    9: [0.0, 0.0, 0.0, 1 / 3, 1 / 3, 1 / 3],         # 6-node triangle
    3: [0.25] * 4,                                   # 4-node quad
    16: [-1 / 12] * 4 + [1 / 3] * 4,                 # 8-node serendipity quad
}
_CORNERS = {1: 2, 8: 2, 2: 3, 9: 3, 3: 4, 16: 4}

# Gauss rules on the reference element, per gmsh boundary element type: [((r, s), weight), ...]
_GL3 = ((-sqrt(0.6), 5 / 9), (0.0, 8 / 9), (sqrt(0.6), 5 / 9))
_LINE_RULE = [((r,), w) for r, w in _GL3]
_QUAD_RULE = [((r, s), wr * ws) for r, wr in _GL3 for s, ws in _GL3]
_a1, _b1, _w1 = 0.059715871789770, 0.470142064105115, 0.132394152788506       # 7-point, degree 5 (Radon)
_a2, _b2, _w2 = 0.797426985353087, 0.101286507323456, 0.125939180544827
_TRI_RULE = [((1 / 3, 1 / 3), 0.5 * 0.225)]
_TRI_RULE += [((r, s), 0.5 * _w1) for r, s in ((_b1, _b1), (_a1, _b1), (_b1, _a1))]
_TRI_RULE += [((r, s), 0.5 * _w2) for r, s in ((_b2, _b2), (_a2, _b2), (_b2, _a2))]
_RULES = {1: _LINE_RULE, 8: _LINE_RULE, 2: _TRI_RULE, 9: _TRI_RULE, 3: _QUAD_RULE, 16: _QUAD_RULE}


def _shape(et: int, r: float, s: float = 0.0):
    """Shape functions N and reference derivatives dN (n_nodes x n_dim) of a gmsh boundary element type,
    in gmsh node order (corners first, then edge mid-nodes in edge order)."""
    if et == 1:                                                   # 2-node line
        N = [0.5 * (1 - r), 0.5 * (1 + r)]
        dN = [[-0.5], [0.5]]
    elif et == 8:                                                 # 3-node line: ends, then middle
        N = [0.5 * r * (r - 1), 0.5 * r * (r + 1), 1 - r * r]
        dN = [[r - 0.5], [r + 0.5], [-2 * r]]
    elif et == 2:                                                 # 3-node triangle
        N = [1 - r - s, r, s]
        dN = [[-1, -1], [1, 0], [0, 1]]
    elif et == 9:                                                 # 6-node triangle
        l1, l2, l3 = 1 - r - s, r, s
        N = [l1 * (2 * l1 - 1), l2 * (2 * l2 - 1), l3 * (2 * l3 - 1), 4 * l1 * l2, 4 * l2 * l3, 4 * l3 * l1]
        dN = [[1 - 4 * l1, 1 - 4 * l1], [4 * l2 - 1, 0], [0, 4 * l3 - 1],
              [4 * (l1 - l2), -4 * l2], [4 * l3, 4 * l2], [-4 * l3, 4 * (l1 - l3)]]
    elif et == 3:                                                 # 4-node quad
        N = [0.25 * (1 - r) * (1 - s), 0.25 * (1 + r) * (1 - s), 0.25 * (1 + r) * (1 + s), 0.25 * (1 - r) * (1 + s)]
        dN = [[-0.25 * (1 - s), -0.25 * (1 - r)], [0.25 * (1 - s), -0.25 * (1 + r)],
              [0.25 * (1 + s), 0.25 * (1 + r)], [-0.25 * (1 + s), 0.25 * (1 - r)]]
    elif et == 16:                                                # 8-node serendipity quad
        N, dN = [], []
        for ri, si in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            N.append(0.25 * (1 + r * ri) * (1 + s * si) * (r * ri + s * si - 1))
            dN.append([0.25 * ri * (1 + s * si) * (2 * r * ri + s * si), 0.25 * si * (1 + r * ri) * (2 * s * si + r * ri)])
        N += [0.5 * (1 - r * r) * (1 - s), 0.5 * (1 + r) * (1 - s * s), 0.5 * (1 - r * r) * (1 + s), 0.5 * (1 - r) * (1 - s * s)]
        dN += [[-r * (1 - s), -0.5 * (1 - r * r)], [0.5 * (1 - s * s), -s * (1 + r)],
               [-r * (1 + s), 0.5 * (1 - r * r)], [-0.5 * (1 - s * s), -s * (1 - r)]]
    else:
        raise ValueError(f"no shape functions for gmsh element type {et}")
    return np.array(N, dtype=float), np.array(dN, dtype=float)


@contextlib.contextmanager
def session(name: str = "model", verbose: bool = False, order: int = 2):
    """Initialise Gmsh with the options used throughout the portfolio, and always finalise."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
        gmsh.option.setNumber("General.NumThreads", 1)      # single-threaded: reproducible meshes
        gmsh.option.setNumber("Mesh.ElementOrder", order)
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)   # 20-node hex, 8-node quad
        gmsh.model.add(name)
        yield gmsh
    finally:
        gmsh.finalize()


def _face_lookup(elements):
    lut = {}
    for etype, elems in elements.items():
        faces = CCX_FACES[etype]
        for eid, nn in elems.items():
            for k, f in enumerate(faces, start=1):
                lut[frozenset(nn[i] for i in f)] = (eid, k)
    return lut


@dataclass
class Mesh:
    nodes: dict[int, tuple[float, float, float]]
    elements: dict[str, dict[int, list[int]]]
    nsets: dict[str, list[int]] = field(default_factory=dict)
    elsets: dict[str, list[int]] = field(default_factory=dict)
    faces: dict[str, list[tuple[int, int]]] = field(default_factory=dict)      # surface -> [(eid, face#)]
    boundary: dict[str, list[tuple[int, list[int]]]] = field(default_factory=dict)  # surface -> [(gmsh type, nodes)]

    # ---- bookkeeping
    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def n_elements(self) -> int:
        return sum(len(e) for e in self.elements.values())

    def n_dofs(self, per_node: int = 3) -> int:
        return per_node * self.n_nodes

    def summary(self) -> str:
        types = ", ".join(f"{k}:{len(v)}" for k, v in self.elements.items())
        return f"{self.n_nodes} nodes, {self.n_elements} elements [{types}], sets {list(self.nsets)}"

    def with_element_type(self, mapping: dict[str, str]) -> "Mesh":
        """Same connectivity, different CalculiX type names (e.g. {'C3D20': 'C3D20R'})."""
        new = {mapping.get(k, k): v for k, v in self.elements.items()}
        return replace(self, elements=new)

    def write_inp(self, path) -> Path:
        return _ccx.write_mesh_inp(path, self.nodes, self.elements, self.nsets, self.elsets)

    def coords(self, ids=None) -> np.ndarray:
        ids = sorted(self.nodes) if ids is None else ids
        return np.array([self.nodes[int(i)] for i in ids])

    def select_nodes(self, predicate) -> list[int]:
        """Node ids for which predicate(x, y, z) is true."""
        return [n for n, (x, y, z) in self.nodes.items() if predicate(x, y, z)]

    def bounding_box(self):
        c = self.coords()
        return c.min(axis=0), c.max(axis=0)

    # ---- loads and boundary helpers
    def _measure(self, et: int, row: list[int]) -> float:
        p = np.array([self.nodes[n] for n in row[:_CORNERS[et]]], dtype=float)
        if len(p) == 2:
            return float(np.linalg.norm(p[1] - p[0]))
        area = 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
        if len(p) == 4:
            area += 0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[0]))
        return float(area)

    def surface_measure(self, surface: str, thickness: float = 1.0) -> float:
        return sum(self._measure(et, row) for et, row in self.boundary[surface]) * thickness

    def surface_nodal_integral(self, surface: str, values_by_node=None) -> dict[int, float]:
        """Consistent nodal distribution of a surface integral: node -> integral of f N_i over a physical surface,
        using the boundary elements' own shape functions and Gauss quadrature, so curved quadratic faces are
        measured exactly (``surface_measure`` uses corner nodes only). ``values_by_node`` maps node id -> scalar
        f; omit it to distribute the true area. With f = h (T - T_inf) the values are the nodal film loads a
        heat-transfer solver applies; summing them gives the integral itself (``surface_integral``)."""
        out: dict[int, float] = {}
        for et, row in self.boundary[surface]:
            xyz = np.array([self.nodes[n] for n in row], dtype=float)
            f = np.ones(len(row)) if values_by_node is None else np.array([values_by_node[n] for n in row], dtype=float)
            acc = np.zeros(len(row))
            for pt, w in _RULES[et]:
                N, dN = _shape(et, *pt)
                tangents = dN.T @ xyz                                   # 1 or 2 tangent vectors
                dA = np.linalg.norm(tangents[0]) if len(tangents) == 1 else np.linalg.norm(np.cross(tangents[0], tangents[1]))
                acc += w * dA * float(N @ f) * N
            for n, v in zip(row, acc):
                out[n] = out.get(n, 0.0) + float(v)
        return out

    def surface_integral(self, surface: str, values_by_node=None) -> float:
        """Integral of a nodal field over a physical surface (the true area when the field is omitted); see
        ``surface_nodal_integral`` for the quadrature."""
        return float(sum(self.surface_nodal_integral(surface, values_by_node).values()))

    def traction_loads(self, surface: str, traction, thickness: float = 1.0) -> dict[int, np.ndarray]:
        """Consistent nodal forces equivalent to a uniform traction vector on a physical surface."""
        t = np.asarray(traction, dtype=float)
        loads: dict[int, np.ndarray] = {}
        for et, row in self.boundary[surface]:
            f_total = t * self._measure(et, row) * thickness
            for nid, frac in zip(row, _LUMP[et]):
                loads[nid] = loads.get(nid, np.zeros(3)) + frac * f_total
        return loads

    def dload_lines(self, surface: str, pressure: float) -> str:
        return "*DLOAD\n" + "".join(f"{eid}, P{k}, {pressure:.10g}\n" for eid, k in self.faces[surface])

    def film_lines(self, surface: str, t_inf: float, h: float) -> str:
        return "*FILM\n" + "".join(f"{eid}, F{k}, {t_inf:.10g}, {h:.10g}\n" for eid, k in self.faces[surface])

    def dflux_lines(self, surface: str, flux: float) -> str:
        return "*DFLUX\n" + "".join(f"{eid}, S{k}, {flux:.10g}\n" for eid, k in self.faces[surface])

    def surface_lines(self, name: str, surface: str) -> str:
        return f"*SURFACE, NAME={name}, TYPE=ELEMENT\n" + "".join(f"{eid}, S{k}\n" for eid, k in self.faces[surface])


def from_gmsh(dim: int = 3) -> Mesh:
    """Convert the current Gmsh model mesh into a :class:`Mesh` (call after ``gmsh.model.mesh.generate``)."""
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = {int(t): (float(coords[3 * i]), float(coords[3 * i + 1]), float(coords[3 * i + 2]))
             for i, t in enumerate(node_tags)}
    elements: dict[str, dict[int, list[int]]] = {}
    etypes, etags, enodes = gmsh.model.mesh.getElements(dim)
    for et, tags, nn in zip(etypes, etags, enodes):
        if int(et) not in GMSH_TO_CCX:
            raise ValueError(f"unsupported gmsh element type {et}: {gmsh.model.mesh.getElementProperties(et)[0]}")
        name, perm = GMSH_TO_CCX[int(et)]
        arr = np.asarray(nn, dtype=int).reshape(len(tags), -1)[:, perm]
        bucket = elements.setdefault(name, {})
        for t, row in zip(tags, arr):
            bucket[int(t)] = row.tolist()
    lut = _face_lookup(elements)
    nsets, elsets, faces, boundary = {}, {}, {}, {}
    for pdim, ptag in gmsh.model.getPhysicalGroups():
        name = gmsh.model.getPhysicalName(pdim, ptag)
        ntags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(pdim, ptag)
        nsets[name] = sorted(int(t) for t in ntags)
        ents = gmsh.model.getEntitiesForPhysicalGroup(pdim, ptag)
        if pdim == dim:
            ids = []
            for e in ents:
                _, tg, _ = gmsh.model.mesh.getElements(pdim, e)
                for arr in tg:
                    ids.extend(int(t) for t in arr)
            elsets[name] = sorted(ids)
        elif pdim == dim - 1:
            flist, blist = [], []
            for e in ents:
                bt, btags, bnodes = gmsh.model.mesh.getElements(pdim, e)
                for et, tags, nn in zip(bt, btags, bnodes):
                    ncorner = _CORNERS[int(et)]
                    arr = np.asarray(nn, dtype=int).reshape(len(tags), -1)
                    for row in arr:
                        hit = lut.get(frozenset(row[:ncorner].tolist()))
                        if hit is None:
                            continue
                        flist.append(hit)
                        blist.append((int(et), row.tolist()))
            faces[name] = flist
            boundary[name] = blist
    return Mesh(nodes, elements, nsets, elsets, faces, boundary)


# ----------------------------------------------------------------------------- primitive meshers
_FACE_NAMES = ("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX")


def _curve_axis(tag: int) -> int:
    """0/1/2 for a straight curve aligned with x/y/z (uses end points, immune to OCC bbox padding)."""
    pts = gmsh.model.getBoundary([(1, tag)], oriented=False)
    p = [np.array(gmsh.model.getValue(0, t, [])) for _, t in pts]
    return int(np.argmax(np.abs(p[1] - p[0])))


def _tag_box_faces(origin, size, names=_FACE_NAMES):
    """Name the six faces of an axis-aligned box by the position of their centre of mass."""
    tol = 1e-6 * max(size)
    for dim, tag in gmsh.model.getEntities(2):
        c = gmsh.model.occ.getCenterOfMass(dim, tag)
        for axis in range(3):
            if abs(c[axis] - origin[axis]) < tol:
                gmsh.model.addPhysicalGroup(2, [tag], name=names[2 * axis])
                break
            if abs(c[axis] - origin[axis] - size[axis]) < tol:
                gmsh.model.addPhysicalGroup(2, [tag], name=names[2 * axis + 1])
                break


def box_hex(L, W, H, nx, ny, nz, order: int = 2, origin=(0.0, 0.0, 0.0), volume_name="BODY") -> Mesh:
    """Structured hexahedral box (transfinite + recombine) with physical faces XMIN..ZMAX."""
    occ = gmsh.model.occ
    v = occ.addBox(*origin, L, W, H)
    occ.synchronize()
    for dim, tag in gmsh.model.getEntities(1):
        n = (nx, ny, nz)[_curve_axis(tag)]
        gmsh.model.mesh.setTransfiniteCurve(tag, n + 1)
    for dim, tag in gmsh.model.getEntities(2):
        gmsh.model.mesh.setTransfiniteSurface(tag)
        gmsh.model.mesh.setRecombine(2, tag)
    gmsh.model.mesh.setTransfiniteVolume(v)
    _tag_box_faces(origin, (L, W, H))
    gmsh.model.addPhysicalGroup(3, [v], name=volume_name)
    gmsh.option.setNumber("Mesh.ElementOrder", order)
    gmsh.model.mesh.generate(3)
    return from_gmsh(3)


def box_tet(L, W, H, size, order: int = 2, origin=(0.0, 0.0, 0.0), volume_name="BODY") -> Mesh:
    """Unstructured tetrahedral box with physical faces XMIN..ZMAX."""
    occ = gmsh.model.occ
    v = occ.addBox(*origin, L, W, H)
    occ.synchronize()
    _tag_box_faces(origin, (L, W, H))
    gmsh.model.addPhysicalGroup(3, [v], name=volume_name)
    gmsh.option.setNumber("Mesh.MeshSizeMin", size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", size)
    gmsh.option.setNumber("Mesh.ElementOrder", order)
    gmsh.model.mesh.generate(3)
    return from_gmsh(3)


def cylinder_tet(L, D, size, order: int = 2, face_names=("XMIN", "XMAX", "LATERAL"), volume_name="BODY",
                 axis_line: bool = False) -> Mesh:
    """Unstructured tetrahedral cylinder of length L along +x and diameter D, one end disc at x = 0.

    Physical surfaces are named after the position of their centre of mass, in ``face_names`` order:
    the x = 0 disc, the x = L disc, the curved face. With ``axis_line`` a straight curve is embedded on
    the axis so the mesh carries nodes at y = z = 0 (convenient for centre-line profiles)."""
    occ = gmsh.model.occ
    v = occ.addCylinder(0.0, 0.0, 0.0, L, 0.0, 0.0, D / 2.0)
    if axis_line:
        p0, p1 = occ.addPoint(0.0, 0.0, 0.0), occ.addPoint(L, 0.0, 0.0)
        line = occ.addLine(p0, p1)
        occ.fragment([(3, v)], [(1, line)])        # imprints the end points into the end discs
    occ.synchronize()
    tol = 1e-6 * max(L, D)
    for dim, tag in gmsh.model.getEntities(2):
        cx = occ.getCenterOfMass(dim, tag)[0]
        name = face_names[0] if abs(cx) < tol else face_names[1] if abs(cx - L) < tol else face_names[2]
        gmsh.model.addPhysicalGroup(2, [tag], name=name)
    gmsh.model.addPhysicalGroup(3, [v], name=volume_name)
    if axis_line:
        gmsh.model.mesh.embed(1, [line], 3, v)
    gmsh.option.setNumber("Mesh.MeshSizeMin", size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", size)
    gmsh.option.setNumber("Mesh.ElementOrder", order)
    gmsh.model.mesh.generate(3)
    return from_gmsh(3)
