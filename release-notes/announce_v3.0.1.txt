# PyLith v3.0.1 now available

I am pleased to announce the release of PyLith 3.0.1, a finite-element code designed to solve dynamic elastic problems and quasistatic viscoelastic problems in tectonic deformation.

You can download the source code and binaries from
  <https://geodynamics.org/resources/pylith>

Documentation
  <https://pylith.readthedocs.org/en/v3.0.1>

## Release Notes

This is a bug fix release with no new features or changes to the user interface.

* Bug fixes
  * Fix lots of small bugs related to running in parallel
  * Fix several discrepancies among the code, examples, and manual
* Examples
  * Added `examples/subduction-3d` steps 1-4 (included in the manual)
  * Added `examples/troubleshooting-2d` (included in the PyLith v3.0 tutorials but not yet added to the manual)
* Documentation
  * Added instructions for how to remove Apple quarantine attributes
  * Fix LaTeX build of documentation (now available at https://pylith.readthedocs.io)
  * Improved instructions on how to run ParaView Python scripts when starting ParaView from a shortcut
  * Added notes indicating steps of examples are not yet updated for v3.0
  * Fix lots of typos

## Binary packages

* Updated PyLith Parameter Viewer (v2.0.0) for Python 3.

## Known issues

* The default PETSc options provide a computationally expensive preconditioner when solving incompressible elasticity problems in parallel. We expect to have a more optimal preconditioner in the next release.
* You may still encounter a few bugs when running in parallel; they appear to cases with specific partitioning of the mesh in relation to one or more faults.

## Contributors

* Brad Aagaard
* Matthew Knepley
* Charles Williams
