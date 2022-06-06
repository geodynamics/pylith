# New in PyLith Version 3.0.0

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

See [Release Notes](../../intro/release-notes.md) for a summary of features and bug fixes for each release.

% End of file
