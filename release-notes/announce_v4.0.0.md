# PyLith v4.0.0 now available

I am pleased to announce the release of PyLith 4.0.0, a finite-element code designed to solve dynamic elastic problems and quasi-static viscoelastic problems in tectonic deformation.

You can download the source code and binaries from
  <https://geodynamics.org/resources/pylith>

Documentation
  <https://pylith.readthedocs.org/en/v4.0.0>

## Release Notes

### Changes in user parameters

* Changed name of fault Lagrange multiplier field for solution component in Python from `lagrange_fault` to `lagrange_multiplier_fault` to match name of solution field in C++ and output.
* Removed support for importing meshes from LaGriT.

### Other changes

* Change in fault tractions are now included in the fault `data_fields` for prescribed slip.
* Fault and boundary orientation directions are now included in the `info_fields` for simulation output.
* State variables are now included in the default `data_fields` for simulation output.
* The default solver settings use the PETSc proper orthogonal decomposition (POD) methodology for initial guess of solutions to improve convergence.
* Add demonstration of `pylith_powerlaw_gendb` in Step 8 of `examples/reverse-2d`.
* Add demonstration of using poroelasticity with porosity as a state variable to `examples/magma-2d`.
* Switched from CppUnit to Catch2 for the C++ testing framework.
* Improve integration with VSCode for testing and debugging (see Developer Guide).
* Bug fixes
  * Fix errors in KinSrcTimeHistory.py
  * Fix creation of PETSc label for edges when importing Gmsh files. This fixes creation of faults with buried edges for 3D meshes imported from Gmsh.
  * Add containers for solution fields for poroelasticity with faults.
* Update PETSc to version 3.20.2.
* Update Python requirement to version 3.8 or later.
* Update Pyre requirement to version 1.1.0 or later.
* Update SpatialData requirement to version 3.1.0 or later.

### Binary packages

* Update to Python 3.10.10, MPICH 4.1.1, OpenMPI 4.1.5, OpenSL 3.1.0, libffi 3.4.4, curl 8.0.1, CMake 3.26.2, PCRE 10.42, SWIG 4.1.1, Catch2 3.3.2, Sqlite 3410200, libtiff 4.5.0, Proj 9.2.0, HDF5 1.14.0, and NetCDF 4.9.2.

### Known issues

* The default PETSc options provide a computationally expensive preconditioner when solving incompressible elasticity problems in parallel. We expect to have a more optimal preconditioner in the next release.

## Contributors

* Brad Aagaard
* Matthew Knepley
* Charles Williams
* Grant Block
* Daniel Douglas
* Lorraine Hwang
* Rezgar Shakeri
* Robert Walker
