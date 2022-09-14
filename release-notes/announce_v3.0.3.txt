# PyLith v3.0.3 now available

I am pleased to announce the release of PyLith 3.0.3, a finite-element code designed to solve dynamic elastic problems and quasistatic viscoelastic problems in tectonic deformation.

You can download the source code and binaries from
  <https://geodynamics.org/resources/pylith>

Documentation
  <https://pylith.readthedocs.org/en/v3.0.3>

## Release Notes

This is a bug fix release with no new features or changes to the user interface.
This list includes changes in both v3.0.2 and v3.0.3.

* Add check of PyLith version against version requirements specified in metadata of parameter files.
* Update defaults to better match most use cases.
  * Use nonlinear solver.
  * Basis order is 1 for solution fields.
  * Basis order is 0 for Cauchy stress and strain.
  * Use ML algebraic multigrid preconditioner (from Trilinos) instead of GAMG preconditioner for more robust solves. This is a temporary change until we find better GAMG settings.
* Update PETSc to v3.17.3.
* Remove obsolete LaTeX documentation.
* Bug fixes
  * Add `viz` directory missing from `examples/subduction-2d` in source distribution.
  * Project output fields using correct PETSc routine (`DMProjectFieldLabel()`). Fixes memory access bugs in both serial and parallel.
  * Fix build warnings.
  * Fix reordering that causes errors when importing Gmsh files.
  * Creating constraints on buried fault edges failed for some mesh distribution cases.
  * Green's function problems did not manage fault impulses on multiple processes.
  * Creating a point mesh for `OutputSolnPoints` failed when running in parallel.
  * PetscSF inconsistencies generated errors at various times when running in parallel.
* Documentation
  * Add discussion of translating boundary value problem information to parameter settings. Add more code blocks to manual.
  * Add discussion of `examples/troubleshooting-2d` to manual.

### Binary packages

* Added PyQT5 Python module for interactive plotting with matplotlib.
* Update PyLith Parameter Viewer to v2.0.1 (fix errors in packaging).
* Update to Python 3.10.6.
* Use `gmforker` process manager with MPICH to avoid localhost name issues.

### Known issues

* The default PETSc options provide a computationally expensive preconditioner when solving incompressible elasticity problems in parallel. We expect to have a more optimal preconditioner in the next release.

## Contributors

* Brad Aagaard
* Matthew Knepley
* Charles Williams
