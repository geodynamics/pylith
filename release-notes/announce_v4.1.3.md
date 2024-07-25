# PyLith v4.1.3 now available

I am pleased to announce the release of PyLith 4.1.3, a finite-element code designed to solve dynamic elastic problems and quasi-static viscoelastic problems in tectonic deformation.

You can download the source code and binaries from
  <https://geodynamics.org/resources/pylith>

Documentation
  <https://pylith.readthedocs.org/en/v4.1.3>

## Release Notes

* **Added**
  * Expanded the section of the manual discussing mesh generation with Cubit and Gmsh. Added a list of useful functions in the Cubit and Gmsh Python interfaces.
  * Added a test to verify Gmsh files with entities that are not embedded generate errors when read via `MeshIOPetsc`.
  * Expanded section of manual discussing installation on Windows Subsystem for Linux, including components needed to run GUI applications (Gmsh and `pylith_viz`).
  * Add Gmsh Python script and corresponding mesh in `examples/magma-2d`.
* **Fixed**
  * Added output of fault traction change when computing static Green's functions.
  * Turn on PETSc HDF5 collective input and output by default. This resolves errors when using HDF5 1.14 on NFS filesystems.
  * Include all necessary files from `examples/subduction-3d` in source and binary packages.
  * Fix typos in documentation and update Components section of User Guide (sync with code).

## Version 4.1.2 (2024/06/12)

* **Fixed**
  * Fix inconsistency in normal direction on fault surfaces. Orientation was correct but direction was flipped at some locations. This affected local slip direction and the resulting deformation close to the fault.

## Version 4.1.1 (2024/06/09)

* **Fixed**
  * Improved setup of variable block Jacobian preconditioner used for elasticity with faults to reduce runtime.
  * Fix several typos and update a few diagrams in the documentation.

### Version 4.1.0 (2024/06/06)

* **Changed**
  * Improved the default preconditioners for poroelasticity for improved scalability.
  * Improved the default preconditioners for elasticity with faults for improved scalability.
  * Replaced use of Cubit journal files with Cubit Python scripts for several examples, and use the skeleton sizing function to increase cell size with distance from the fault.
  * Removed ParaView Python scripts (replaced by `pylith_viz`)
* **Added**
  * New 2D and 3D crustal strike-slip faults examples based on the 2019 Ridgecrest earthquake.
  * New 2D subduction zone outer-rise faulting example examining poroelastic response to bending stresses.
  * New `pylith_viz` utility for plotting PyLith results.
  * Updated `examples/strikeslip-2d` and `examples/reverse-2d` to demonstrate use of uniform refinement and higher order discretizations to improve resolution of solution.
  * Documentation
    * Developer Guide: Added description of the process for adding event logging and evaluating performance, including performance logging.
    * Developer Guide: Added checklist for what is needed when contributing examples.
* **Fixed**
  * Updated `examples/subduction-3d` Steps 6, 7, and 8 to work with PyLith v3+.
  * Fix performance bottleneck when reading Gmsh files.
  * Remove extra (wrong) kernels for poroelasticity when using state variables.
  * Update Python unit tests setup for compatibility with the current unittest API (use `load_tests()`).

### Known issues

* The new default preconditioner for simulations using elasticity and faults can cause convergence issues when running in parallel in which fault faces lie on the boundaries between processors. The workaround is to use the previous preconditioner provided in `share/settings/solver_fault_fieldsplit.cfg`.
* The default PETSc options provide a computationally expensive preconditioner when solving incompressible elasticity problems. We expect to have a more optimal preconditioner in the next release.

## Contributors

* Brad Aagaard ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0000-0002-8795-9833](https://orcid.org/0000-0002-8795-9833)
* Matthew Knepley ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0000-0002-2292-0735](https://orcid.org/0000-0002-2292-0735)
* Charles Williams ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0000-0001-7435-9196](https://orcid.org/0000-0001-7435-9196)
* Daniel Douglas ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0000-0002-7871-018X](https://orcid.org/0000-0002-7871-018X)
* Evan Marschall ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0009-0003-0916-6656](https://orcid.org/0009-0003-0916-6656)
* Zechao Zhuo ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0000-0002-8163-5731](https://orcid.org/0000-0002-8163-5731)
