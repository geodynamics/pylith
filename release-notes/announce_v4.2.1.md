# PyLith v4.2.1 now available

I am pleased to announce the release of PyLith 4.2.1, a finite-element code designed to solve dynamic elastic problems and quasi-static viscoelastic problems in tectonic deformation.

You can download the source code and binaries from
  <https://geodynamics.org/resources/pylith>

Documentation
  <https://pylith.readthedocs.org/en/v4.2.1>

## Release Notes

* **Changed**
  * Update PETSc to v3.23.5
  * Use the VTU (XML) format for VTK files instead of the legacy ASCII format (required for PETSc v3.23.5).
* **Added**
  * Improved documentation for `pylith_eqinfo` and illustrate use in `examples/strikdslip-2d`, `examples/crustal-strikeslip-2d`, and `examples/crustal-strikeslip-3d`.
  * Improved the organization of the governing equations section and added documentation for poroelasticity with prescribed fault slip, including governing equations and default PETSc solver settings. Also added a full-scale test.
* **Fixed**
  * Account for processes without cells in initial condition patch when verifying configuration.
  * Avoid divide by zero for `KinSrcRamp` when final slip is zero.
  * Remove `solid_bulk_modulus` as a spatial database field for poroelasticity; compute it from the other fields.
  * Fixed typos in setting up gravity with poroelasticity.
  * Update default solver settings for poroelasticity and faults.
  * Fix pythia import in `pylith_eqinfo`.

### Known issues

* The new default preconditioner for simulations using elasticity and faults can cause convergence issues when running in parallel in which fault faces lie on the boundaries between processors. The workaround is to use the previous preconditioner provided in `share/settings/solver_fault_fieldsplit.cfg`.
* The default PETSc options provide a computationally expensive preconditioner when solving incompressible elasticity problems. We expect to have a more optimal preconditioner in the next release.

## Contributors

* Brad Aagaard [ORCID iD 0000-0002-8795-9833](https://orcid.org/0000-0002-8795-9833)
* Matthew Knepley [ORCID iD 0000-0002-2292-0735](https://orcid.org/0000-0002-2292-0735)
* Charles Williams [ORCID iD 0000-0001-7435-9196](https://orcid.org/0000-0001-7435-9196)
