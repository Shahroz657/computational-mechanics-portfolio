"""fealib.plotting - one consistent, print-friendly matplotlib style for every figure in the portfolio."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from cycler import cycler  # noqa: E402

# categorical palette (fixed order, colour-blind checked): blue, orange, aqua, yellow, magenta, green, violet, red
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def setup() -> None:
    plt.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 170, "font.size": 10, "font.family": "sans-serif",
        "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK2,
        "axes.labelcolor": INK, "axes.titlecolor": INK, "axes.titleweight": "bold", "axes.titlesize": 11,
        "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "lines.linewidth": 2.0, "lines.markersize": 6.5, "legend.frameon": False, "legend.fontsize": 9,
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.prop_cycle": cycler(color=PALETTE),
    })


def save(fig, path, **kw) -> None:
    fig.savefig(path, bbox_inches="tight", **kw)
    plt.close(fig)


# ----------------------------------------------------------------------------- 2D field plots
def triangulation_2d(mesh):
    """matplotlib Triangulation of a 2D fealib.mesh.Mesh (corner nodes only) + the node-id order used."""
    import matplotlib.tri as mtri
    import numpy as np

    ids = sorted(mesh.nodes)
    idx = {n: i for i, n in enumerate(ids)}
    xy = np.array([mesh.nodes[n][:2] for n in ids])
    tris = []
    for etype, elems in mesh.elements.items():
        ncorner = 3 if etype[-1] in "36" else 4
        for nn in elems.values():
            c = [idx[k] for k in nn[:ncorner]]
            tris.append(c[:3])
            if ncorner == 4:
                tris.append([c[0], c[2], c[3]])
    return mtri.Triangulation(xy[:, 0], xy[:, 1], np.array(tris)), ids


def sequential_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUE)


def field_2d(ax, mesh, values_by_node: dict, *, levels=24, cmap=None, mesh_lines=False, label=None):
    """Filled contour of a nodal field on a 2D mesh. ``values_by_node`` maps node id -> scalar."""
    import numpy as np

    tri, ids = triangulation_2d(mesh)
    vals = np.array([values_by_node[n] for n in ids])
    cs = ax.tricontourf(tri, vals, levels=levels, cmap=cmap or sequential_cmap())
    if mesh_lines:
        element_outlines(ax, mesh)
    ax.set_aspect("equal")
    ax.grid(False)
    cb = ax.figure.colorbar(cs, ax=ax, shrink=0.85, pad=0.02)
    if label:
        cb.set_label(label)
    return cs


def element_outlines(ax, mesh, color=INK2, lw=0.3, alpha=0.45):
    """Draw the true element edges of a 2D mesh (quads stay quads)."""
    from matplotlib.collections import PolyCollection

    polys = []
    for etype, elems in mesh.elements.items():
        ncorner = 3 if etype[-1] in "36" else 4
        for nn in elems.values():
            polys.append([mesh.nodes[k][:2] for k in nn[:ncorner]])
    ax.add_collection(PolyCollection(polys, facecolor="none", edgecolor=color, linewidth=lw, alpha=alpha))
