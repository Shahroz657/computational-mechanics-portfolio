# Expert review of the AI-generated FEA solution

**Verdict: reject.** The submission satisfies its own validation step (reaction force = applied load) and
still reports a maximum stress that is 43 % below the hand calculation, a deflection that is two million
times too small, and a safety factor built on the wrong strength property. Every point below is
reproduced by [`run.py`](run.py); the numbers are from the committed run
([`results/summary.json`](results/summary.json)).

## What the submission claims

From [`results/ai_report.md`](results/ai_report.md): 794 linear tetrahedra (10 mm), clamped face fixed,
500 N on one tip node; maximum von Mises stress **85.8 MPa** at node 10 (x = 0, y = 10, z = 20 mm);
tip deflection **0.0000 mm ("negligible")**; safety factor **5.01** against S_ut = 430 MPa; "the reaction
force equals the applied load, confirming the model is correct."

## Reference values (closed form + converged FE)

| Quantity | Value | Source |
|---|---|---|
| Root bending stress `σ = M c / I = P L (h/2) / (b h³/12)` | **150.0 MPa** | hand calculation |
| Bending stress at x = 2H = 40 mm | 135.0 MPa | hand calculation; converged C3D20R model gives 135.0 MPa |
| Tip deflection (Timoshenko, κ = 0.8497) | **4.008 mm** | hand calculation; converged C3D20R model gives 3.992 mm |
| Static safety factor, ductile steel | **1.83** = S_y / σ_root = 275 / 150 | yield criterion |

## Findings

| # | Severity | Finding | Evidence |
|---|---|---|---|
| 1 | **Blocker** | **Inconsistent units.** Geometry in mm and E = 200e9 (Pa). In an mm–N system E must be 200 000 MPa; the stiffness is therefore 1000× too high. Stresses in a statically determinate linear problem do not depend on E, so this error is invisible in the stress output and shows only in displacement: 1.7e-6 mm reported vs 4.008 mm expected. The "negligible deflection" conclusion is an artefact of the unit error (× ≈ 1000) compounded by over-stiff elements (× ≈ 2). | check C1 |
| 2 | **Blocker** | **The reported maximum is not converged and is read at a singular location.** The maximum sits on the clamped face, where a fully fixed boundary meets a free surface: the stress there grows without limit under refinement. Refining from 10 mm to 2.5 mm moves the reported maximum from 85.8 to 138.5 MPa (+62 %), with no sign of levelling off. The submission calls a single-mesh result "confirmed". | checks C3, C4; figure below |
| 3 | **Blocker** | **Wrong strength criterion.** For a ductile steel under static load the safety factor is taken against yield (275 MPa), not ultimate (430 MPa). Combined with the under-predicted stress, the reported margin of 5.01 should be 1.83. That is the difference between "comfortable" and "marginal for a load that is probably not known to ± 10 %". | check C7 |
| 4 | Major | **Linear tetrahedra in bending.** C3D4 elements are far too stiff in bending: on the 10 mm mesh the top-fibre stress at x = 2H is 32.9 MPa where theory gives 135 MPa, and even 49 556 elements (2.5 mm) reach only 119.7 MPa. Project 01 shows quadratic elements hit 135.0 MPa with 20 elements. | check C6, refinement table |
| 5 | Major | **Point load on a solid mesh.** A single nodal force on a continuum is a stress singularity by construction; the von Mises stress at that node rises from 13.6 to 56.5 MPa across the refinement series and would rise indefinitely. The resultant belongs on the end face as a traction (consistent nodal loads) or through a distributing coupling. | check C8 |
| 6 | Major | **Validation logic.** Equilibrium of reactions is necessary, not sufficient: it holds exactly for every mesh in the refinement series, including the ones that are 60 % wrong. A convergence check and a comparison with `M c / I` were the minimum required, and both were available for free. | check C2 (passes) vs C3, C5 |
| 7 | Minor | **No stated element type, integration or solver settings** beyond "tetrahedral"; the mesh count is given but not the element order. A reviewer cannot judge the model without this. | report text |

## What "reliable" would have looked like

The corrected model in `run.py` uses consistent mm–N–MPa units, 1 280 C3D20R hexahedra with four
elements through the depth, the 500 N resultant applied as a uniform traction on the end face, and
reports the stress **away from the clamp singularity**: 135.0 MPa at x = 2H against 135.0 MPa from
`M c / I`, tip deflection 3.992 mm against 4.008 mm (Timoshenko), and a safety factor of 1.83 on yield.
Project 01 documents the full convergence study behind that choice of element and mesh.

![singularity](results/singularity.png)

## Objective checks used (executable, `review_checks()` in `run.py`)

| Check | Result |
|---|---|
| C1 unit consistency: E / σ_max within 1e2 … 1e5 for a metal at working stress | FAIL (2.3e9) |
| C2 global equilibrium: reactions equal applied load within 1e-6 | PASS |
| C3 mesh convergence of the reported quantity: < 5 % change under 4× refinement | FAIL (+62 %) |
| C4 reported maximum not at a point load or constraint edge | FAIL (clamp face) |
| C5 hand-calculation plausibility: σ_max within 25 % of `M c / I`, deflection within 20 % of theory | FAIL |
| C6 element formulation suited to bending (no linear tets) | FAIL |
| C7 static safety factor of a ductile part based on yield | FAIL |
| C8 load applied as traction or coupling, not a single nodal force | FAIL |

Seven of eight checks fail. Checks C1–C5 are generic and can be run on any submitted static analysis
without knowing the answer in advance; C6–C8 are modelling-practice checks that need an engineer to
read the deck.
