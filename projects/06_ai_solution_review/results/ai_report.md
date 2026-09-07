# FEA result: cantilever bar under 500 N tip load

**Model.** Linear static analysis in CalculiX. The bar (400 x 20 x 20 mm) was meshed with
794 tetrahedral elements (10 mm). The clamped face was fully fixed
(all DOF) and the 500 N load was applied at the tip node. Steel: E = 200 GPa, nu = 0.3.

**Results.**
- Maximum von Mises stress: **85.8 MPa** at node 10 (x = 0 mm, y = 10 mm, z = 20 mm).
- Tip deflection: **0.0000 mm** (negligible).
- Safety factor against ultimate strength (430 MPa): **5.01**.

**Validation.** The reaction force at the clamp is 500.0 N, which equals the applied
load, confirming the model is correct. The design is safe with a comfortable margin.
