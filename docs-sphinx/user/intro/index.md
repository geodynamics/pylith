# Introduction

## Overview

PyLith is portable, scalable software for simulation of crustal deformation across spatial scales ranging from meters to hundreds of kilometers and temporal scales ranging from milliseconds to thousands of years.
Its primary applications are quasistatic and dynamic modeling of earthquake faulting.

## New in PyLith Version 3.0.0dev

* Major rewrite of the finite-element implementation to support higher order discretizations and flexible specification of the governing equations.
    * Use of pointwise functions to implement governing equations;
    * Higher order discretizations;
    * Problem specification independent of cell shape (quadrilateral vs triangle, hexahedron vs tetrahedron);
    * Incompressible elasticity;
    * Poroelasticity;
    * Use of PETSc time-stepping algorithms;
* Simulations now require metadata, such as description, command line arguments, and PyLith version compatibility.
* New utilities
    * `pyre_doc.py` Display facilities and components available for a Pyre component.
    * `pylith_cfgsearch` Find files matching criteria for metadata.
    * `pylith_runner` Run all simulations in a specified path.
* New examples
    * Simple 2-D and 3-D examples of Dirichlet and Neumann boundary conditions without faults;
    * Prescribed slip on a 2-D through-going strike-slip fault;
    * Gravitational body forces with elasticity and incompressible elasticity;
    * Distributed surface loads using Neumann boundary conditions; and
    * Prescribed slip on a reverse fault with a splay fault;
* Developer guide has moved to <https://pylith.readthedocs.io>.
* Updated to Python 3
    * Pythia/Pyre, spatialdata, and PyLith have all been migrated to Python 3.
    * The nemesis package has been merged into Pyre/Pyre.

See [Release Notes](release-notes.md) for a summary of features and bug fixes for each release.

## History

PyLith 1.0 was the first version to allow the solution of both implicit
(quasistatic) and explicit (dynamic) problems and was a complete rewrite of
the original PyLith (version 0.8). PyLith 1.0 combines the functionality of
EqSim {cite}`Aagaard:etal:2001a`, {cite}`Aagaard:etal:2001b` and PyLith 0.8.
PyLith 0.8 was a direct descendant of LithoMop and was the first version that
ran in parallel, as well as providing several other improvements over
LithoMop. LithoMop was the product of major re-engineering of Tecton, a
finite-element code for simulating static and quasistatic crustal deformation.
The major new features present in LithoMop included dynamic memory allocation
and the use of the Pyre simulation framework and PETSc solvers. EqSim was
written by Brad Aagaard to solve problems in earthquake dynamics, including
rupture propagation and seismic wave propagation.

The release of PyLith 1.0 has been followed by additional releases that expand
the number of features as well as improve performance. The PyLith 1.x series
of releases allows the solution of both quasistatic and dynamic problems in
one, two, or three dimensions. The code runs in either serial or parallel, and
the design allows for relatively easy scripting using the Python programming
language. Material properties and values for boundary and fault conditions are
specified using spatial databases, which permit easy prescription of complex
spatial variations of properties and parameters. Simulation parameters are
generally specified through the use of simple ASCII files or the command line.
At present, mesh information may be provided using a simple ASCII file (PyLith
mesh ASCII format) or imported from CUBIT or LaGriT, two widely-used meshing
packages. The elements currently available include a linear bar in 1D, linear
triangles and quadrilaterals in 2D, and linear tetrahedra and hexahedra in 3D.
Materials presently available include isotropic elastic, linear Maxwell
viscoelastic, generalized Maxwell viscoelastic, power-law viscoelastic, and
Drucker-Prager elastoplastic. Boundary conditions include Dirichlet
(prescribed displacements and velocities), Neumann (traction), point forces,
and absorbing boundaries. Cohesive elements are used to implement slip across
interior surfaces (faults) with both kinematically-specified fault slip and
slip governed by fault constitutive models. PyLith also includes an interface
for computing static Green's functions for fault slip.

PyLith 2.0 replaced the finite-element data structures provided by the C++
Sieve implementation with those provided by the C DMPlex implementation. The
newly developed DMPlex implementation by the PETSc developers conforms to the
PETSc data manager (DM) interface, thereby providing tighter integration with
other PETSc data structures, such as vectors and matrices. Other improvements
include significantly reduced memory use and memory balancing.

PyLith 3.0 involves restructuring the code to permit more flexible
specifications of governing equations and discretization and use of PETSc
time-stepping algorithms. The finite-element integrations, constraints, and
transformations are done through a suite of point-wise functions. The PETSc
DMPlex interface calls these functions to perform the finite-element
integrations.

PyLith is under active development and we expect a number of additions and
improvements in the near future. Likely enhancements will include additional
bulk and fault constitutive models, coupled quasistatic and dynamic
simulations for earthquake cycle modeling, and coupling between elasticity,
heat flow, and/or fluid flow.

## PyLith Workflow

PyLith is one component in the process of investigating problems in tectonics
({numref}`fig:workflow:summary`). Given a geological problem of
interest, a scientist must first provide a geometrical representation of the
desired structure. Once the structure has been defined, a computational mesh
must be created. PyLith presently provides three mesh importing options: CUBIT
Exodus format, LaGriT GMV and Pset files, and PyLith mesh ASCII format. The
modeling of the physical processes of interest is performed by a code such as
PyLith. Present output consists of VTK or HDF5/Xdmf files which can be used by
a number of visualization codes (e.g., ParaView, Visit, and Matlab).

:::{figure-md} fig:workflow:summary
<img src="figs/workflow.*" alt="Workflow involved in going from geologic structure to problem analysis." width="100%" />

Workflow involved in going from geologic structure to problem analysis.
:::
