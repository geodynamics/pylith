# PyLith v3.0.0 now available

I am pleased to announce the release of PyLith 3.0.0, a finite-element code designed to solve dynamic elastic problems and quasistatic viscoelastic problems in tectonic deformation.

This major release is the result of more than 5 years of work to overhaul the code in support of a wide range of new features.

You can download the source code and binaries from
  <https://geodynamics.org/resources/pylith>

Documentation
  <https://pylith.readthedocs.org/en/v3.0.0>

## Migrating from version 2.X TO 3.0.0

Version 3.0.0 includes major changes to the underlying finite-element formulation and implementation in order to support a more flexible specification of the governing equations and higher order basis functions. These changes affect how simulations are defined. Parameter files for previous versions will need to be updated; the changes are too complex for a simple translation table. Some features present in v2.2.2, such as spontaneous rupture and finite strain, have not yet been implemented in the new formulation.

## Release Notes

* Multiphysics
  * Elasticity for linear isotropic materials and linear Maxwell, generalized Maxwell, and power law viscoelastic models
  * Incompressible elasticity for linear isotropic materials
  * Prescribed slip for quasistatic and dynamic simulations
* Higher order basis functions
    Allow user to select order of basis functions independent of the mesh (which defines the geometry). This permits higher resolution for a given mesh.
* Switch to using PETSc time-stepping (TS) algorithms
  Replace simple Python-based time-stepping implementations with PETSc time-stepping algorithms that provide support for higher order discretization in time and real adaptive time stepping.
* Static Green's functions with user-specified discretization of fault slip impulses
* Import finite-element meshes from Cubit (Exodus II), Gmsh, and LaGriT
* Modular approach for initial conditions
* Output of subfields with user-defined basis order
* Simulation metadata with command line utility for searching metadata
* Convert to Python 3
* Convert LaTeX documentation to Sphinx + MyST
* Testing with the Method of Manufactured Solutions
* Automatically assign label value for fault cohesive cells (`id` setting is obsolete).
* Use `description` for descriptive labels and `label` and `label_value` for tagging entities. PyLith's use of`label` and `label_value` now corresponds to PETSc labels and label values.

### Deprecated features

We plan to discontinue support for reading LaGriT mesh files in version 3.2. Gmsh provides an open-source alternative with a graphical user interface.

## Known issues

* Running in parallel has a few minor bugs due to communication mismatches and over-aggressive error checking. We will be fixing these in the next week.
* We will be updating the 3D subduction zone example (examples/subduction-3d) to v3.0.0 in the next week, including providing the input mesh file; in the meantime do not attempt to run this example.
* We have included Gmsh in the binary packages.
For Linux there are additional libraries that must be installed for Gmsh to run; these are associated with the graphical user interface and included in most default installations.

## Contributors

* Brad Aagaard
* Matthew Knepley
* Charles Williams
* Robert Walker
* Chris Mills
* Shengduo Liu
* Thea Ragon
* Alex Berne
* Jed Brown
* Rey Koki
* Kali Allison
* Lorraine Hwang
