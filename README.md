# Computational mechanics portfolio · Muhammad Shahroz

[![verify](https://github.com/Shahroz657/computational-mechanics-portfolio/actions/workflows/verify.yml/badge.svg)](https://github.com/Shahroz657/computational-mechanics-portfolio/actions/workflows/verify.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![site](https://img.shields.io/badge/portfolio-site-1f5fb4.svg)](https://shahroz657.github.io/computational-mechanics-portfolio/)

Scripted finite-element engineering in Python: geometry from code (CadQuery, Gmsh OCC), meshes from
the Gmsh API, plain-text CalculiX decks, solver runs from the command line, results parsed into NumPy,
and **every model verified against a closed-form solution with executable PASS/FAIL checks**. No
graphical interface is used anywhere in the pipeline. The quick projects re-run in GitHub Actions on
every push.

Portfolio site with figures and a summary of each study: **https://shahroz657.github.io/computational-mechanics-portfolio/**

## Projects

| # | Project | Physics | Toolchain | Verified result (from the committed run) |
|---|---|---|---|---|
| 01 | [Cantilever verification and element convergence](projects/01_cantilever_verification) | linear static, 3D solids | Gmsh structured hex/tet → CalculiX | 22 runs; converged C3D20R tip deflection 0.9964 × Timoshenko, bending stress 135.0 vs 135.0 MPa; shear locking (C3D8) and hourglassing (C3D8R, 27.7 × too flexible) reproduced |
| 02 | [Plate with a hole: Kt sweep](projects/02_plate_with_hole_stress_concentration) | 2D plane stress, stress concentration | Gmsh size fields → CalculiX CPS8 | Kt within 1.4–2.0 % of Peterson's chart for d/W = 0.1–0.6, peak at θ = 0°, converged to 0.09 % |
| 03 | [Modal and buckling analysis](projects/03_modal_and_buckling) | vibration, linear buckling | CalculiX `*FREQUENCY`, `*BUCKLE` | 12 runs; 10 modes classified automatically, first weak-axis mode 51.06 vs 50.96 Hz, every bending mode within 0.31 % of exact Timoshenko; first buckling load 5150 vs 5140 N (+0.19 %); a CalculiX `*BUCKLE` conditioning pitfall isolated and documented |
| 04 | [Pin-fin heat transfer](projects/04_pin_fin_heat_transfer) | steady conduction + convection | Gmsh OCC cylinder → CalculiX heat transfer | 10 solves; centre-line temperature within 0.014 °C of Incropera's convective-tip solution, heat rate within 0.22 % as printed and 0.015 % after correcting CalculiX's reaction-flux bookkeeping (traced node by node), energy balance closed to 9e-5 % |
| 05 | [L-bracket design optimisation](projects/05_bracket_design_optimization) | static strength, stiffness, mass | CadQuery → STEP → Gmsh → CalculiX → SciPy | 41 solves: 35-point DOE, response surfaces (R² 0.999), SLSQP optimum t = 7.19 mm, R = 10 mm verified by FE at SF 2.56 and 0.590 mm, 26.9 % lighter than the thick-plate baseline; fillet stress converged to 0.8 % while the fixed-edge singularity grew 21 % |
| 06 | [Reviewing an AI-generated FEA solution](projects/06_ai_solution_review) | review, singularities, units | CalculiX, executable checks | 7 of 8 objective checks fail on the AI solution: units (deflection 2·10⁶ × too small), non-converged singular maximum (+62 % under refinement), wrong strength criterion (SF 5.01 → 1.83) |
| 07 | [Fatigue life estimation](projects/07_fatigue_life_estimation) | stress-life fatigue | NumPy, own ASTM E1049 rainflow | rainflow count identical to the reference package in every bin (1082.5 cycles/block); four mean-stress corrections span a factor of 5.05 in life; constant-amplitude chain reproduces N(S) exactly |

Each project folder has a `run.py` that regenerates everything in `results/` (figures, CSV tables,
`summary.json` with a `checks` dictionary) and exits non-zero if any check fails, plus a README that
states the model, the reference equations, the numbers and the limitations.

## How every study is built

```
CadQuery / Gmsh OCC  ->  Gmsh API mesh  ->  CalculiX deck (text)  ->  ccx (CLI)  ->  .frd / .dat parsers
        geometry         2nd-order hex/tet     BCs, loads, steps       solver         NumPy arrays
                                                                                          |
                     closed-form reference  <-------------------------------------  verify + checks
```

The shared package [`fealib/`](fealib) is small on purpose and fully readable:

- [`fealib/mesh.py`](fealib/mesh.py) – Gmsh session helper, structured/unstructured box meshers, conversion of
  any Gmsh model to CalculiX element types with the node-order permutations (C3D20, C3D10, C3D15, CPS8…),
  physical groups → node/element sets, element-face lists for pressure/film loads, **consistent nodal loads**
  for uniform tractions on quadratic faces
- [`fealib/ccx.py`](fealib/ccx.py) – deck text helpers, solver runner (fails loudly on `*ERROR`), `.frd`
  parser (any nodal block: DISP, STRESS, NDTEMP…), `.dat` parser (eigenfrequencies, buckling factors,
  reaction totals), von Mises / principal stresses
- [`fealib/analytic.py`](fealib/analytic.py) – Timoshenko/Euler beams, cantilever modes, Euler buckling,
  Peterson and Kirsch stress concentration, Incropera fins, Lamé cylinders
- [`tests/`](tests) – unit tests for the parsers, the mesh conversion and the closed forms, plus an
  end-to-end CalculiX run compared with Timoshenko

## Running it yourself

```bash
git clone https://github.com/Shahroz657/computational-mechanics-portfolio
cd computational-mechanics-portfolio
python -m venv .venv && source .venv/bin/activate      # Python 3.10-3.12
pip install -e ".[dev]"                                # numpy, scipy, gmsh, cadquery, matplotlib, pytest ...
python -m pytest -q
cd projects/01_cantilever_verification && python run.py
```

CalculiX (`ccx`) is found on `PATH`, in `$CCX`, or inside a FreeCAD bundle:

| Platform | Command |
|---|---|
| Ubuntu / Debian | `sudo apt install calculix-ccx libglu1-mesa libxrender1 libxcursor1 libxft2 libxinerama1` |
| macOS | `brew install calculix-ccx` or `export CCX=/Applications/FreeCAD.app/Contents/Resources/bin/ccx` |
| any | set `CCX=/path/to/ccx` |

Gmsh comes from the `gmsh` wheel on PyPI (no separate install). Versions used for the committed
results: CalculiX 2.22, Gmsh 4.13/4.15, CadQuery 2.8, Python 3.12, macOS arm64.

## Earlier CAD work

Before this repository my design work lived in SolidWorks and PTC Creo: the Ventus-I agricultural
hexacopter (design lead, seven-person team, IMechE UAS Challenge), a V6 twin-turbo engine assembly
(Creo, team of four), and personal modelling and rendering projects. Renders are on the
[portfolio site](https://shahroz657.github.io/computational-mechanics-portfolio/#cad).

## About

B.Sc. Mechanical Engineering, National University of Sciences and Technology (NUST), 2024, with a thesis on
environmental ageing of hybrid composites (ASTM D7264 / D2344, 72 specimens, 1 000+ measurements).
MIT Emerging Talent certificate in computer and data science (2025). Since 2024 I design and review
engineering evaluation tasks for AI-training data teams (micro1, Mercor, Turing), and I am extending the
thesis with a computational model of hygrothermal ageing in epoxy laminates ahead of PhD applications.

[LinkedIn](https://www.linkedin.com/in/shahroz657) · [GitHub](https://github.com/Shahroz657) · muhammadshahroz1019@gmail.com

## License

MIT. Result figures and tables in `projects/*/results` are produced by the scripts and may be reused
with attribution.
