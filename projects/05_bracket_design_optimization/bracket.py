"""Parametric L-bracket built with CadQuery (millimetres).

The geometry is code, so every design point of the optimisation is regenerated from two numbers:
plate thickness ``t`` and inner fillet radius ``R``. The model is exported to STEP for meshing and
its exact volume is used for the mass objective.
"""
from __future__ import annotations

from pathlib import Path

import cadquery as cq

WB = 50.0            # bracket width (y)
HV = 80.0            # height of the vertical (mounting) plate (z)
LA = 100.0           # arm length beyond the plate face (x)
HOLE_D = 9.0         # two mounting holes through the plate
HOLE_Z = HV - 20.0


def bracket(t: float, R: float) -> cq.Workplane:
    plate = cq.Workplane("XY").box(t, WB, HV, centered=False)          # x: 0..t, y: 0..WB, z: 0..HV
    arm = cq.Workplane("XY").box(t + LA, WB, t, centered=False)        # x: 0..t+LA, z: 0..t
    body = plate.union(arm)
    body = body.edges(cq.selectors.NearestToPointSelector((t, WB / 2, t))).fillet(R)   # inner corner
    holes = (cq.Workplane("YZ").pushPoints([(WB / 4, HOLE_Z), (3 * WB / 4, HOLE_Z)])
             .circle(HOLE_D / 2).extrude(t))
    return body.cut(holes)


def export_step(t: float, R: float, path) -> Path:
    path = Path(path)
    cq.exporters.export(bracket(t, R), str(path))
    return path


def export_svg(t: float, R: float, path, view=(-1.0, -1.6, 0.9)) -> Path:
    path = Path(path)
    cq.exporters.export(bracket(t, R), str(path), opt={"projectionDir": view, "width": 640, "height": 420,
                                                       "showAxes": False, "strokeWidth": 0.6, "showHidden": False})
    return path


def volume_mm3(t: float, R: float) -> float:
    return float(bracket(t, R).val().Volume())


def mass_kg(t: float, R: float, rho_kg_m3: float = 7850.0) -> float:
    return volume_mm3(t, R) * 1e-9 * rho_kg_m3
