Conductivity Solver v1.2 (27.4.26)
A numerical solver for the calculation of thermal and electrical conductivity in heterogenous media.
Full tensor support and arbitrary geometry handling.

Features
1. Full Thermal Conductivity Tensor
* **Periodic Boundary Conditions (PBC):** Implementation of PBCs in all spatial directions.
* **Full Tensor Support:** Unlike previous versions limited to Dirichlet/Neumann conditions, v1.2 can now compute the **complete effective thermal conductivity tensor** $\lambda_{eff}$.
* Improved accuracy for anisotropic rock samples and complex pore networks.

2. Arbitrary Geometries
* **Beyond Cubic Domains:** Support for non-cuboid geometries (e.g., spheres, irregular rock cuttings).
* **Backward and Forward Search:** A newly implemented routine ensures that the solver correctly identifies domain boundaries in non-standard geometries.
* *Reference:* This methodology is based on the research presented in [Insert Your Paper Name Here].

3. Solver Optimization & CSR Migration
To accommodate the complexity of periodic boundaries and larger datasets, the internal matrix logic has been refactored:
* **CSR Matrix Format:** Transitioned from diagonal vector storage to **Compressed Sparse Row (CSR)** format for the system matrix.
* **Refactored Linear Algebra:** All core functions for matrix-vector multiplication and algebra have been optimized for the CSR format.

----------------------------------------------------------------------------------------------------------------------------------
Version 1.1 (24.01.24):
  - code is organized in seperate files 
  - added flux field computation for postprocessing 
  - added automated paraview load file for postprocessing 
  - adjustments for user-friendliness

  
- preprocessing tutorial video: https://youtu.be/RUQssngar3Y
- postprocessing tutorial video: https://youtu.be/coOlP0cAj5g
