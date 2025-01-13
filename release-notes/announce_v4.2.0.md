# PyLith v4.2.0 now available

I am pleased to announce the release of PyLith 4.2.0, a finite-element code designed to solve dynamic elastic problems and quasi-static viscoelastic problems in tectonic deformation.

You can download the source code and binaries from
  <https://geodynamics.org/resources/pylith>

Documentation
  <https://pylith.readthedocs.org/en/v4.2.0>

## Release Notes

* **Changed**
  * Updated `examples/strikeslip-2d` Steps 4-7 to use a more realistic slip distribution and mesh refinement in output.
  * Updated to PETSc 3.22.0
  * Switch CI from Azure Pipelines to GitHub Actions.
* **Added**
  * Default filenames for progress monitor and parameters file are set from the simulation name like
  the other output files.
  * Consistency check for material properties and scales used in nondimensionalization (currently just the shear modulus).
  * Added section to User Guide on troubleshooting solver issues.
  * Added section to User Guide on how to run PyLith binary offline.
  * Allow output on a finer resolution mesh than used in the simulation to facilitate accurate representation of fields with a basis order of 2 or greater.
* **Fixed**
  * Fixed inconsistency in normal direction on fault surfaces. Orientation was correct but direction was flipped at some locations. This affected local slip direction and the resulting deformation close to the fault. This bug fix was not in version 4.1.3.
  * Updated autoconf numpy macros for compatibility with location of include files in numpy version 2.x.

### Known issues

* The new default preconditioner for simulations using elasticity and faults can cause convergence issues when running in parallel in which fault faces lie on the boundaries between processors. The workaround is to use the previous preconditioner provided in `share/settings/solver_fault_fieldsplit.cfg`.
* The default PETSc options provide a computationally expensive preconditioner when solving incompressible elasticity problems. We expect to have a more optimal preconditioner in the next release.

## Contributors

* Brad Aagaard ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0000-0002-8795-9833](https://orcid.org/0000-0002-8795-9833)
* Matthew Knepley ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0000-0002-2292-0735](https://orcid.org/0000-0002-2292-0735)
* Charles Williams ![ORCID iD](/_static/images/ORCIDiD_icon32x32.png){w=16px}[0000-0001-7435-9196](https://orcid.org/0000-0001-7435-9196)
