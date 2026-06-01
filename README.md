* [Changelog and Features](#Changelog-and-Features)
* [Tutorial](#Tutorial)
## Tutorial
## Mini-Tutorial: Setting Up a Case

**Prerequisite:** You will need a **Fortran compiler** (e.g., `gfortran`) to run this tool.

Follow these simple steps to set up and run your simulation:

1. **Download Files:** Download all files from this repository and place them together in a single directory.
2. **Configure Settings:** Open the `settings.f90` file in a text editor, adjust the parameters for your specific case, and save the file.
3. **Compile:** Open your terminal in the case directory and compile the code using the following command:
```bash
   gfortran -fopenmp -O3 conductivitySolver.f90 -o solverRun
```
4. **Set Number of CPUs:** Define the number of CPU cores you want to use for the parallel computation (e.g. nCPU=8):
```bash
   export OMP_NUM_THREADS=8
```
5. **Run:**  Execute the compiled program to start the simulation:
```bash
   ./solverRun
```

## Changelog and Features
**Conductivity Solver v1.2 (27.4.26)**

A numerical solver for the calculation of thermal and electrical conductivity in heterogenous media.
Full tensor support and arbitrary geometry handling.

Features
1. Full Thermal Conductivity Tensor
* Periodic Boundary Conditions (PBC): Implementation of PBCs in all spatial directions.
* Full Tensor Support: Unlike previous versions limited to Dirichlet/Neumann conditions, v1.2 can now compute the complete effective conductivity tensor $\lambda_{eff}$.

2. Arbitrary Geometries
* Beyond Cubic Domains: Support for non-cuboid geometries (e.g., spheres, irregular rock cuttings).
* Backward and Forward Search: A newly implemented routine ensures that the solver correctly identifies domain boundaries in non-standard geometries.
* *Reference:* This methodology is based on the research presented in [xxxxxxx].

3. Solver Optimization & CSR Migration: To accommodate the complexity of periodic boundaries and larger datasets, the internal matrix logic has been refactored
* CSR Matrix Format: Transitioned from diagonal vector storage to Compressed Sparse Row (CSR) format for the system matrix.
* Refactored Linear Algebra: All core functions for matrix-vector multiplication and algebra have been optimized for the CSR format.

----------------------------------------------------------------------------------------------------------------------------------
**Conductivity Solver v1.1 (24.01.24):**
  - code is organized in seperate files 
  - added flux field computation for postprocessing 
  - added automated paraview load file for postprocessing 
  - adjustments for user-friendliness

  
- preprocessing tutorial video: https://youtu.be/RUQssngar3Y
- postprocessing tutorial video: https://youtu.be/coOlP0cAj5g
