# PyLith

## Description

PyLith is a finite-element code with a primary focus on modeling interseismic and coseismic deformation of Earth's crust and upper mantle.
PyLith supports 2D and 3D static, quasistatic (neglecting inertia), and dynamic (including inertia) formulations of the governing equations, which can be elasticity, incompressble elasticity, or poroelasticity.
A variety of elastic and viscoelastic bulk rheologies are supported.
Boundary conditions include Dirichlet, Neumann, and absorbing boundaries.
Faults are treated as interior interfaces.
Currently, only the slip must be prescribed (kinematic rupture).
We plan to reimplement spontaneous rupture (fault friction) in an upcoming release; see `ref`{sec-development-plan} for more information.

The code is written in C++ and Python and uses MPI for parallel processing.
We use the Pyre Python framework to setup the simulation.
We leverage [PETSc](https://petsc.org) for finite-element data structures and operations as well as linear and nonlinear solvers.

## Other sources of documentation

* This manual in other formats [epub](https://pylith.readthedocs.io/_/downloads/en/latest/epub/) [pdf](https://pylith.readthedocs.io/_/downloads/en/latest/pdf/)
* [PyLith tutorials](https://geodynamics.org/courses/PyLith)
* [SpatialData documentation](https://spatialdata.readthedocs.io)
* [PyLith installer documentation](https://pylith-installer.readthedocs.io) (for installing from source)

## PyLith development team

The primary developers are:

* Brad Aagaard (U.S. Geological Survey)
* Matthew Knepley (University at Buffalo)
* Charles Williams (GNS Science)

```{include} ../LICENSE.md
```

```{toctree}
---
caption: Table of Contents
hidden: True
---
intro/index.md
user/index.md
developer/index.md
references.md
```

