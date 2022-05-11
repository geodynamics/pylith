# PyLith

```{admonition} Under construction
:class: warning

We are in the process of migrating the documentation from LaTeX producing a PDF file to Markedly Structured Text (MyST) and Sphinx producing this online documentation, a PDF file, and epub.

Most of the user guide has been migrated, but there are still some examples that are missing documentation.
```

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

* PyLith v3.0.0rc1 manual [PDF](https://github.com/geodynamics/pylith/releases/download/v3.0.0beta1/pylith-3.0.0b1_manual.pdf)
* [PyLith tutorials](https://wiki.geodynamics.org/software:pylith:start)
* [PyLith installer documentation](https://pylith_installer.readthedocs.io)
* [Spatial Data documentation](https://spatialdata.readthedocs.io)

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

