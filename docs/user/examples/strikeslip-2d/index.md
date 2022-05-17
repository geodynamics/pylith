# Horizontal Cross-Section of Strike-Slip Fault (2D)

The files are in the directory `examples/strikeslip-2d`.

## Overview

This suite of examples demonstrates some basic concepts of using PyLith to solve the static and quasistatic boundary elasticity equation for a horizontal cross-section of a strike-slip fault ({numref}`fig:example:strikeslip:2d:overview`) with nonuniform material properties.
The fault extends the entire length of the domain.
The shear modulus is larger on the +x side of the fault.
This example builds on the previous examples and adds complexity through a series of steps:

:Step 1: Static coseismic slip with Dirichlet (displacement) boundary conditions.
:Step 2: Quasistatic coseismic slip with time-dependent Dirichlet (displacement) boundary conditions.
:Step 3: Quasistatic slip with two ruptures and time-dependent Dirichlet (displacement) boundary conditions.

:::{figure-md} fig:example:strikeslip:2d:overview
<img src="figs/geometry.*" alt="Diagram of geometry for strike-slip fault." scale="75%"/>

Diagram of geometry for domain with a strike-slip fault.
The domain extends from -50 km to +50 km in the x direction and from -75 km to +75 km in the y direction.
We refer to the domain boundaries using the names shown in the diagram.
:::

:::{important}
We decribe how to generate the finite-element mesh using both Gmsh and Cubit.
The files for both methods are included.
For Step 1 we provide PyLith parameter files for both meshes; for Steps 2 and 3 we only provide the Parameter files that use the Gmsh file.
:::

## Example Workflow

:::{toctree}
meshing-gmsh.md
meshing-cubit.md
common-information.md
step01-slip.md
step02-slip-velbc.md
step03-multislip-velbc.md
exercises.md
:::
