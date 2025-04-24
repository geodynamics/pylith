# Significant Recent Changes

See [Release Notes](../../intro/release-notes.md) for a summary of features and bug fixes for each release.

## Version 5.0

* Use the VTU (XML) format for VTK files instead of the legacy ASCII format.
* Renamed `VertexGroup` to `BoundaryGroup` in `meshio.gmsh_utils` and changed default behavior to not be recursive (generate "face" groups, not "vertex" groups).
* Material `description` property is no longer used; a deprecation warning is printed to stdout if it is specified. This feature will be removed in v6.0.

## Version 4.2

* Default filenames for progress monitor and parameters file are set from the simulation name like the other output files.
* Allow output on a finer resolution mesh than used in the simulation to facilitate accurate representation of fields with a basis order of 2 or greater.
* Fixed inconsistency in normal direction on fault surfaces. Orientation was correct but direction was flipped at some locations. This affected local slip direction and the resulting deformation close to the fault. This bug fix in a PETSc branch was not in version 4.1.3.

## Version 4.1

* Improved runtime performance, including better preconditioners for elasticity with faults and poroelasticity
* 2D and 3D crustal strike-slip faults examples based on the 2019 Ridgecrest earthquake.
* 2D subduction zone outer-rise faulting example examining poroelastic response to bending stresses.
* `pylith_viz` utility for plotting PyLith results.
* Updated `examples/strikeslip-2d` and `examples/reverse-2d` to demonstrate use of uniform refinement and higher order discretizations to improve resolution of solution.

## Version 4.0

* Changed name of fault Lagrange multiplier field for solution component in Python from `lagrange_fault` to `lagrange_multiplier_fault` to match name of solution field in C++.
* Removed support for importing meshes from LaGriT.
* Change in fault tractions are now included in the fault `data_fields` for prescribed slip.
* Fault and boundary orientation directions are now included in the `info_fields` for simulation output.
* State variables are now included in the default `data_fields` for simulation output.
* The default solver settings use the PETSc proper orthogonal decomposition (POD) methodology for initial guess of solutions to improve convergence.

## Version 3.0

* Major rewrite of the finite-element implementation to support higher order discretizations and flexible specification of the governing equations.
  * Use of pointwise functions to implement governing equations;
  * Higher order discretizations;
  * Problem specification independent of cell shape (quadrilateral vs triangle, hexahedron vs tetrahedron);
  * Incompressible elasticity;
  * Poroelasticity; and
  * Use of PETSc time-stepping algorithms.
* Simulations now require metadata, such as description, command line arguments, and PyLith version compatibility.
* New utilities
  * `pyre_doc.py` Display facilities and components available for a Pyre component;
  * `pylith_cfgsearch` Find files matching criteria for metadata; and
  * `pylith_runner` Run all simulations in a specified path.
* New examples
  * Simple 2-D and 3-D examples of Dirichlet and Neumann boundary conditions without faults;
  * Prescribed slip on a 2-D through-going strike-slip fault;
  * Gravitational body forces with elasticity and incompressible elasticity;
  * Distributed surface loads using Neumann boundary conditions; and
  * Prescribed slip on a reverse fault with a splay fault.
* Documentation is now available online at <https://pylith.readthedocs.io>.
* Import finite-element meshes from Gmsh in addition to Cubit (Exodus II) and LaGriT.
* Updated to Python 3.
  * Pythia/Pyre, spatialdata, and PyLith have all been migrated to Python 3; and
  * The nemesis package has been merged into Pyre/Pyre.
