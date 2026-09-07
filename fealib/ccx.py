"""fealib.ccx - a small, explicit driver for CalculiX CrunchiX (ccx).

* locate the solver (``$CCX``, PATH, FreeCAD bundle, common prefixes)
* write mesh include files (``*NODE`` / ``*ELEMENT`` / ``*NSET`` / ``*ELSET``)
* run the solver and fail loudly on ``*ERROR``
* parse ASCII ``.frd`` result files into numpy arrays (DISP, STRESS, NDTEMP, ...)
* parse ``.dat`` files (eigenfrequencies, buckling factors, set totals)

Decks stay plain text so every model in this repository can be opened and audited.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_CANDIDATES = [
    "ccx", "ccx_2.22", "ccx_2.21", "ccx_2.20",
    "/Applications/FreeCAD.app/Contents/Resources/bin/ccx",
    "/opt/homebrew/bin/ccx", "/usr/local/bin/ccx", "/usr/bin/ccx",
]


def find_ccx() -> str:
    """Return the path of the ccx executable (``$CCX`` wins, then PATH, then known bundles)."""
    env = os.environ.get("CCX")
    if env:
        return env
    for cand in _CANDIDATES:
        hit = shutil.which(cand)
        if hit:
            return hit
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    raise FileNotFoundError(
        "CalculiX solver 'ccx' not found. Set CCX=/path/to/ccx or install it "
        "(macOS: FreeCAD bundle or `brew install calculix-ccx`; Ubuntu: `apt install calculix-ccx`)."
    )


# ----------------------------------------------------------------------------- deck text
def _chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def fmt_nodes(nodes: dict[int, tuple[float, float, float]], nset: str = "NALL") -> str:
    out = [f"*NODE, NSET={nset}"]
    for nid in sorted(nodes):
        x, y, z = nodes[nid]
        out.append(f"{nid}, {x:.10g}, {y:.10g}, {z:.10g}")
    return "\n".join(out) + "\n"


def fmt_elements(eltype: str, elems: dict[int, list[int]], elset: str = "EALL") -> str:
    out = [f"*ELEMENT, TYPE={eltype}, ELSET={elset}"]
    for eid in sorted(elems):
        rows = list(_chunks([eid, *elems[eid]], 10))
        for k, row in enumerate(rows):
            out.append(", ".join(map(str, row)) + ("," if k < len(rows) - 1 else ""))
    return "\n".join(out) + "\n"


def fmt_set(kind: str, name: str, ids) -> str:
    key = "NSET" if kind.upper() == "NSET" else "ELSET"
    out = [f"*{key}, {key}={name}"]
    out.extend(", ".join(map(str, row)) for row in _chunks(sorted({int(i) for i in ids}), 12))
    return "\n".join(out) + "\n"


def fmt_cload(loads: dict[int, np.ndarray], tol: float = 0.0) -> str:
    """``*CLOAD`` block from a {node: [Fx,Fy,Fz]} dictionary (zero components are skipped)."""
    out = ["*CLOAD"]
    for nid in sorted(loads):
        for dof, val in enumerate(np.asarray(loads[nid], dtype=float), start=1):
            if abs(val) > tol:
                out.append(f"{nid}, {dof}, {val:.10g}")
    return "\n".join(out) + "\n"


def fmt_boundary(target, first_dof: int, last_dof: int | None = None, value: float | None = None) -> str:
    """One ``*BOUNDARY`` line for a node set name or an iterable of node ids."""
    last_dof = first_dof if last_dof is None else last_dof
    tail = f", {value:.10g}" if value is not None else ""
    if isinstance(target, str):
        return f"*BOUNDARY\n{target}, {first_dof}, {last_dof}{tail}\n"
    return "*BOUNDARY\n" + "".join(f"{n}, {first_dof}, {last_dof}{tail}\n" for n in sorted(target))


def write_mesh_inp(path, nodes, elements: dict[str, dict[int, list[int]]], nsets=None, elsets=None) -> Path:
    """Write a mesh include file: nodes, one ``*ELEMENT`` block per type (all in ELSET=EALL), then sets."""
    parts = [fmt_nodes(nodes)]
    for eltype, elems in elements.items():
        parts.append(fmt_elements(eltype, elems, elset="EALL"))
        if len(elements) > 1:
            parts.append(fmt_set("ELSET", f"E_{eltype}", elems.keys()))
    for name, ids in (nsets or {}).items():
        parts.append(fmt_set("NSET", name, ids))
    for name, ids in (elsets or {}).items():
        parts.append(fmt_set("ELSET", name, ids))
    path = Path(path)
    path.write_text("".join(parts))
    return path


def material_block(name: str, E: float | None = None, nu: float | None = None, rho: float | None = None,
                   k: float | None = None, cp: float | None = None, alpha: float | None = None) -> str:
    out = [f"*MATERIAL, NAME={name}"]
    if E is not None:
        out += ["*ELASTIC", f"{E:.10g}, {nu:.10g}"]
    if rho is not None:
        out += ["*DENSITY", f"{rho:.10g}"]
    if k is not None:
        out += ["*CONDUCTIVITY", f"{k:.10g}"]
    if cp is not None:
        out += ["*SPECIFIC HEAT", f"{cp:.10g}"]
    if alpha is not None:
        out += ["*EXPANSION", f"{alpha:.10g}"]
    return "\n".join(out) + "\n"


# ----------------------------------------------------------------------------- running
@dataclass
class RunResult:
    job: str
    workdir: Path
    frd: Path
    dat: Path
    log: str
    seconds: float


def run_ccx(inp, ccx: str | None = None, nproc: int | None = None, timeout: float = 3600) -> RunResult:
    """Run ``ccx -i <job>`` in the deck's directory. Raises on non-zero exit or any ``*ERROR`` line."""
    inp = Path(inp).resolve()
    job = inp.stem
    exe = ccx or find_ccx()
    env = dict(os.environ)
    nthreads = str(nproc or os.cpu_count() or 1)
    env.setdefault("OMP_NUM_THREADS", nthreads)
    env.setdefault("CCX_NPROC_STIFFNESS", nthreads)
    env.setdefault("CCX_NPROC_RESULTS", nthreads)
    t0 = time.time()
    proc = subprocess.run([exe, "-i", job], cwd=inp.parent, env=env, capture_output=True, text=True, timeout=timeout)
    log = proc.stdout + proc.stderr
    (inp.parent / f"{job}.log").write_text(log)
    if proc.returncode != 0 or "*ERROR" in log:
        errs = "\n".join(l for l in log.splitlines() if "ERROR" in l.upper()) or log[-3000:]
        raise RuntimeError(f"ccx failed for {job} (exit {proc.returncode}):\n{errs}")
    frd = inp.with_suffix(".frd")
    if not frd.exists():
        raise RuntimeError(f"ccx produced no {frd.name}; log tail:\n{log[-2000:]}")
    return RunResult(job, inp.parent, frd, inp.with_suffix(".dat"), log, time.time() - t0)


# ----------------------------------------------------------------------------- .frd parsing
@dataclass
class ResultBlock:
    """One nodal result block (e.g. DISP, STRESS, NDTEMP) of one step/mode."""
    name: str
    components: list[str]
    node_ids: np.ndarray
    values: np.ndarray          # shape (n_nodes, n_components)
    step_value: float = float("nan")
    step: int = 0
    _index: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        self._index = {int(n): i for i, n in enumerate(self.node_ids)}

    def at(self, nodes) -> np.ndarray:
        """Values for the given node id(s), in the order given."""
        if np.isscalar(nodes):
            return self.values[self._index[int(nodes)]]
        return self.values[[self._index[int(n)] for n in nodes]]

    def column(self, comp: str) -> np.ndarray:
        return self.values[:, self.components.index(comp)]

    def as_dict(self) -> dict[int, np.ndarray]:
        return {int(n): self.values[i] for i, n in enumerate(self.node_ids)}


@dataclass
class FrdResult:
    nodes: dict[int, tuple[float, float, float]]
    blocks: list[ResultBlock]

    def get(self, name: str, step: int | None = None) -> ResultBlock:
        hits = [b for b in self.blocks if b.name == name]
        if not hits:
            raise KeyError(f"no result block {name!r}; have {[b.name for b in self.blocks]}")
        if step is None:
            return hits[-1]
        return hits[step]

    def all(self, name: str) -> list[ResultBlock]:
        return [b for b in self.blocks if b.name == name]

    def coords(self, ids) -> np.ndarray:
        return np.array([self.nodes[int(i)] for i in ids])


def _floats12(line: str, start: int) -> list[float]:
    vals = []
    for k in range(start, len(line) - 11, 12):
        chunk = line[k:k + 12]
        if chunk.strip():
            vals.append(float(chunk))
    return vals


def read_frd(path) -> FrdResult:
    """Parse an ASCII CalculiX ``.frd`` file (nodes + all nodal result blocks, in file order)."""
    lines = Path(path).read_text().splitlines()
    nodes: dict[int, tuple[float, float, float]] = {}
    blocks: list[ResultBlock] = []
    i, n = 0, len(lines)
    step_value, step_no = float("nan"), 0
    while i < n:
        line = lines[i]
        if line.startswith("    2C"):                     # node block
            i += 1
            while i < n and lines[i].startswith(" -1"):
                l = lines[i]
                nodes[int(l[3:13])] = (float(l[13:25]), float(l[25:37]), float(l[37:49]))
                i += 1
            continue
        if line.startswith("    3C"):                     # element block: skip
            i += 1
            while i < n and not lines[i].startswith(" -3"):
                i += 1
            i += 1
            continue
        if line.lstrip().startswith("100C"):              # step header
            k = line.index("100C")
            try:
                step_value = float(line[k + 10:k + 22])
            except ValueError:
                step_value = float("nan")
            try:
                step_no = int(line[k + 57:k + 62])
            except ValueError:
                step_no += 1
            i += 1
            continue
        if line.startswith(" -4"):                        # result block
            name = line[5:13].strip()
            i += 1
            comps = []
            while i < n and lines[i].startswith(" -5"):
                comps.append(lines[i][5:13].strip())
                i += 1
            ids, rows = [], []
            while i < n and (lines[i].startswith(" -1") or lines[i].startswith(" -2")):
                l = lines[i]
                if l.startswith(" -1"):
                    ids.append(int(l[3:13]))
                    rows.append(_floats12(l, 13))
                else:
                    rows[-1].extend(_floats12(l, 13))
                i += 1
            values = np.array(rows, dtype=float)
            if values.ndim == 2:
                comps = comps[:values.shape[1]]
            blocks.append(ResultBlock(name, comps, np.array(ids, dtype=int), values, step_value, step_no))
            continue
        i += 1
    return FrdResult(nodes, blocks)


def von_mises(stress: np.ndarray) -> np.ndarray:
    """Von Mises equivalent stress from rows of [SXX, SYY, SZZ, SXY, SYZ, SZX] (CalculiX order)."""
    s = np.atleast_2d(stress)
    sxx, syy, szz, sxy, syz, szx = (s[:, k] for k in range(6))
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3.0 * (sxy ** 2 + syz ** 2 + szx ** 2))


def principal_stresses(stress: np.ndarray) -> np.ndarray:
    """Principal stresses (descending) for rows of [SXX,SYY,SZZ,SXY,SYZ,SZX]."""
    s = np.atleast_2d(stress)
    out = np.empty((s.shape[0], 3))
    for r, (sxx, syy, szz, sxy, syz, szx) in enumerate(s):
        t = np.array([[sxx, sxy, szx], [sxy, syy, syz], [szx, syz, szz]])
        out[r] = np.sort(np.linalg.eigvalsh(t))[::-1]
    return out


# ----------------------------------------------------------------------------- .dat parsing
_NUM = r"[-+]?\d+\.\d+E[-+]\d+|[-+]?\d+\.?\d*(?:E[-+]?\d+)?"


@dataclass
class DatResult:
    eigenfrequencies_hz: list[float]
    eigenvalues: list[float]
    buckling_factors: list[float]
    totals: list[dict]      # {"label": "force", "set": "FIXED", "time": 1.0, "values": [..]}


def read_dat(path) -> DatResult:
    """Extract eigenfrequencies, buckling factors and ``TOTALS=ONLY`` sums from a ``.dat`` file."""
    text = Path(path).read_text()
    lines = text.splitlines()
    freqs, eigs, bucks, totals = [], [], [], []
    section = None
    for idx, line in enumerate(lines):
        u = line.upper()
        if "E I G E N V A L U E" in u and "O U T P U T" in u:
            section = "eig"
            continue
        if "B U C K L I N G" in u and "O U T P U T" in u:
            section = "buck"
            continue
        if "O U T P U T" in u or "F A C T O R S" in u or "M A S S" in u:
            section = None
        toks = line.split()
        if section == "eig" and len(toks) == 5 and toks[0].isdigit():
            eigs.append(float(toks[1]))
            freqs.append(float(toks[3]))
        elif section == "buck" and len(toks) == 2 and toks[0].isdigit():
            bucks.append(float(toks[1]))
        m = re.match(r"\s*total (\w[\w ]*?) (?:\([^)]*\) )?for set (\S+) and time\s+(\S+)", line)
        if m:
            j = idx + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            vals = [float(v) for v in lines[j].split()] if j < len(lines) else []
            totals.append({"label": m.group(1).strip(), "set": m.group(2), "time": float(m.group(3)), "values": vals})
    return DatResult(freqs, eigs, bucks, totals)
