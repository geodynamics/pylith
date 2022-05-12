(sec-examples-2d-box)=
# Axial and Shear Deformation (2D Box)

The files are in the directory `examples/box-2d`.

## Overview

This suite of examples demonstrates some basic concepts of using PyLith to solve the static and quasistatic boundary elasticity equation in a 2D box ({numref}`fig:example:box:2d:overview`) with uniform material properties.
This example incrementally adds complexity through a series of steps:

:Step 1: Axial extension with Dirichlet (displacement) boundary conditions.
:Step 2: Shear deformation with Dirichlet (displacement) boundary conditions.
:Step 3: Shear deformation with Dirichlet (displacement) and Neumann (traction) boundary conditions.
:Step 4: Same as Step 2 but with initial conditions equal to the analytical solution.
:Step 5: Shear deformation with time-dependent Dirichlet (displacement) and Neumann (traction) boundary conditions.

:::{figure-md} fig:example:box:2d:overview
<img src="figs/geometry.*" alt="Diagram of geometry for 2D box." scale="75%"/>

Diagram of geometry for 2D box that extends from -6 km to +6 km in the x direction and from -16 km to 0 km in the y direction.
We refer to the domain boundaries using the names shown in the diagram.
:::

## Example Workflow

:::{toctree}
meshing.md
common-information.md
step01-axialdisp.md
step02-sheardisp.md
step03-sheardisptract.md
step04-sheardispic.md
step05-sheardisptractrate.md
exercises.md
:::
