(sec:example:subduction:2d)=
# Vertical Cross-Section of Subduction Zone (2D)

## Overview

This example examines quasistatic interseismic and coseismic deformation in 2D for a subduction zone (see {numref}`fig:example:subduction:2d:overview`).
It is based on the 2011 M9.0 Tohoku earthquake off the east coast of Japan.
{numref}`fig:example:subduction:2d:steps` shows the three steps of increasing complexity.
Step 1 focuses on the coseismic slip, Step 2 focuses on interseismic deformation, and Step 3 combines the two into a pseudo-earthquake cycle deformation simulation.
Step 4 focuses on using the change in tractions from Step 1 to construct a simulation with afterslip controlled by frictional sliding.
Steps 5 and 6 replace the prescribed aseismic slip on the subducting slab in Step 2 with a frictional interface, producing spontaneous earthquake ruptures and creep.

:::{figure-md} fig:example:subduction:2d:overview
<img src="figs/cartoon_general.*" alt="Cartoon of subduction zone example." width="100%"/>

Cartoon of subduction zone example.
:::

:::{figure-md} fig:example:subduction:2d:steps
<img src="figs/steps.*" alt="Diagram of fault slip and boundary conditions for each step in the subduction zone example." width="100%" />

Diagram of fault slip and boundary conditions for each step in the subduction zone.
:::

## Features

:::{admonition} TODO @brad
:class: error

Use script to create `features.md`.
:::

PyLith features discussed in this example:

* Static solution
* Quasi-static solution
* CUBIT/Trelis mesh generation w/APREPRO
* Nonplanar geometry
* Variable mesh resolution
*  Linear triangular cells
*  HDF5 output
*  Dirichlet displacement and velocity boundary conditions
*  ZeroDispDB spatial database
*  UniformDB spatial database
*  SimpleDB spatial database
*  SimpleGridDB
*  Multiple materials
*  Nonlinear solver
*  Plane strain linearly elastic material
*  Plane strain linear Maxwell viscoelastic material
*  Prescribed slip
*  Spontaneous rupture
*  Multiple faults
*  Spatially variable coseismic slip
*  Spatially variable aseismic creep
*  Afterslip via fault friction
*  Static friction
*  Slip-weakening friction
*  Rate-state friction

All of the files necessary to run the examples are contained in the directory `examples/2d/subduction`.

## Example Workflow

:::{toctree}
meshing.md
common-information.md
step01-coseismic-slip.md
step02-interseismic-deformation.md
step03-psuedo-earthquake-cycle-model.md
step04-frictional-afterslip.md
step05-spontaneous-slip-weakening.md
step06-spontaneous-rate-state.md
step07-twofaults-maxwell.md
step08-twofaults-powerlaw.md
:::

## Exercises

The list below includes some suggested modifications to these examples that will allow you to become more familiar with PyLith while examining some interesting physics.

*  Change the resolution of the mesh by editing the `mesh_tri3.jou` journal file.
Change the resolution and bias factor.
*  Add depth dependent viscosity to the mantle and crust.
This requires using the linear Maxwell plane strain bulk constitutive model in the crust as well and creating spatial databases that include viscosity for the crust.
Specifying a depth dependent variation in the parameters will require adding points, updating num-locs accordingly, and changing data-dim to 1.
*  Modify the spatial database files for the material properties to use depth-dependent elastic properties based on PREM (Dziewonski and Anderson, 1981, 10.1016/0031-9201(81)90046-7). See <http://ds.iris.edu/ds/products/emc-prem/> for a simple table of values. Add points, update num-locs accordingly, and change data-dim to 1.
*  Modify the CUBIT journal files to use quad4 cells rather than tri3 cells.
This requires using the pave mesh scheme.
*  Modify Steps 5 and 6 to use a user-defined variable time step.
Experiment with longer time steps between earthquake ruptures and smaller time steps around the time of the earthquake ruptures.
Can you develop a simple algorithm for choosing the time step?
*  Adjust the parameters of the friction models and examine the effects on the deformation and the convergence of the nonlinear solve.
In which cases do you need to adjust the time step to retain reasonable convergence?
