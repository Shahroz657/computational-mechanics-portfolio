| quantity                  | ai                              | corrected                                                           |
|:--------------------------|:--------------------------------|:--------------------------------------------------------------------|
| max stress reported (MPa) | 85.8 (at node 10)               | 150.0 nominal at the root (M c / I); FE at 2H 135.0 vs theory 135.0 |
| tip deflection (mm)       | 0.00000                         | 3.992 (Timoshenko 4.008)                                            |
| safety factor             | 5.01 (S_ut / sigma)             | 1.83 (S_y / sigma_root)                                             |
| elements                  | 794 C3D4, single mesh           | 1280 C3D20R, convergence shown in project 01                        |
| load application          | single nodal force              | uniform traction via consistent nodal loads                         |
| units                     | mm geometry with E = 200e9 (Pa) | mm-N-MPa: E = 200 000 MPa                                           |