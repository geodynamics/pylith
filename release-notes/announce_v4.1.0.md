# PyLith v4.1.0 now available

I am pleased to announce the release of PyLith 4.1.0, a finite-element code designed to solve dynamic elastic problems and quasi-static viscoelastic problems in tectonic deformation.

You can download the source code and binaries from
  <https://geodynamics.org/resources/pylith>

Documentation
  <https://pylith.readthedocs.org/en/v4.1.0>

## Release Notes

* **Changed**
  * Improved the default preconditioners for poroelasticity for improved scalability.
  * Improved the default preconditioners for elasticity with faults for improved scalability.
  * Replaced use of Cubit journal files with Cubit Python scripts for several examples, and use the skeleton sizing function to increase cell size with distance from the fault.
  * Removed ParaView Python scripts (replaced by `pylith_viz`)
* **Added**
  * New 2D and 3D crustal strike-slip faults examples based on the 2019 Ridgecrest earthquake.
  * New 2D subduction zone outer-rise faulting example examining poroelastic response to bending stresses.
  * New `pylith_viz` utility for plotting PyLith results.
  * Updated `examples/strikeslip-2d` and `examples/reverse-2d` to demonstrate use of uniform refinement and higher ordder discretizations to improve resolution of solution.
  * Documentation
    * Developer Guide: Added description of the process for adding event logging and evaluating performance, including performance logging.
    * Developer Guide: Added checklist for what is needed when contributing examples.
* **Fixed**
  * Updated `examples/subduction-3d` Steps 6, 7, and 8 to work with PyLith v3+.
  * Fix performance bottleneck when reading Gmsh files.
  * Remove extra (wrong) kernels for poroelasticity when using state variables.
  * Update Python unit tests setup for compatibility with the current unittest API (use `load_tests()`).

### Known issues

* The default PETSc options provide a computationally expensive preconditioner when solving incompressible elasticity problems. We expect to have a more optimal preconditioner in the next release.

## Contributors

* Brad Aagaard ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0000-0002-8795-9833](https://orcid.org/0000-0002-8795-9833)
* Matthew Knepley ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0000-0002-2292-0735](https://orcid.org/0000-0002-2292-0735)
* Charles Williams ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0000-0001-7435-9196](https://orcid.org/0000-0001-7435-9196)
* Daniel Douglas ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0000-0002-7871-018X](https://orcid.org/0000-0002-7871-018X)
* Evan Marschall ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0009-0003-0916-6656](https://orcid.org/0009-0003-0916-6656)
* Zechao Zhuo ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0000-0002-8163-5731](https://orcid.org/0000-0002-8163-5731)
