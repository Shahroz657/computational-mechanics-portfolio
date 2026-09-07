"""fealib - small, explicit helpers for scripted FEA: Gmsh meshing, CalculiX decks, result parsing,
and closed-form reference solutions used to verify every model in this portfolio."""
from . import analytic, ccx, mesh, plotting  # noqa: F401
__all__ = ["analytic", "ccx", "mesh", "plotting"]
