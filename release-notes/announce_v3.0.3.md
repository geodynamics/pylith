# PyLith v3.0.3 now available

I am pleased to announce the release of PyLith 3.0.3, a finite-element code designed to solve dynamic elastic problems and quasi-static viscoelastic problems in tectonic deformation.

You can download the source code and binaries from
  <https://geodynamics.org/resources/pylith>

Documentation
  <https://pylith.readthedocs.org/en/v3.0.3>

## Release Notes

This is a bug fix release with no new features or changes to the user interface.

* Fixed duplicate integration of fault terms if a fault had one material on one side and multiple materials on the other side.
* Fixed bugs related to running in parallel.
  * Creating constraints on buried fault edges failed for some mesh distribution cases.
  * Green's function problems did not manage fault impulses on multiple processes.
  * Creating a point mesh for `OutputSolnPoints` failed when running in parallel.
  * PetscSF inconsistencies generated errors at various times when running in parallel.
* Update to PETSc 3.18.0.

**Note**: We now use PETSc routines to write the HDF5 files.
As a result, there is one change to the layout: `topology/cells` is now `viz/topology/cells`.
The corresponding Xdmf files reflect this change.

### Binary packages

* Update to Python 3.10.6.
* Use `gmforker` process manager with MPICH to avoid localhost name issues.

### Known issues

* The default PETSc options provide a computationally expensive preconditioner when solving incompressible elasticity problems in parallel. We expect to have a more optimal preconditioner in the next release.

## Contributors

* Brad Aagaard
* Matthew Knepley
* Charles Williams
