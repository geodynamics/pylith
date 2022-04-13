# Vertical Cross-Section of a Reverse Fault with Splay (2D)

## Overview

This suite of examples demonstrates use of a number of features for a simple two-dimensional model.
This example also shows how to produce a mesh with a somewhat complex geometry.
Although the problem geometry ({numref}`fig:example:2d:reverse:geometry`) includes a simple planar splay fault intersecting a planar thrust fault, the first 3 steps actually focus on gravitational body forces, reference stresses, and incompressible elasticity.
The fourth example demonstrates the use of traction boundary conditions to represent a surface load.
The remainder of the examples focus on slip on one or more faults, including an example of multiple ruptures on a single fault.
To keep the meshing and computation time in these examples short, we limit our model to a 200 km x 100 km domain and we will use a relatively coarse discretization.

:::{figure-md} fig:example:2d:reverse:geometry
<img src="figs/geometry.*" alt="Geometry used for 2D reverse fault example." width="100%"/>

Geometry used for 2D reverse fault example.
:::

Note that although we label the different parts of the mesh as slab, crust, and wedge, the actual thrust fault only extends 60 km downdip, and we do not provide bottom boundaries for the crust and slab.
The files associated with this suite of examples are contained in the directory `examples/2d/reverse`.
This directory contains several files:

* **\*.jou** Files used to construct the finite-element mesh using CUBIT/Trelis.
* **\*.spatialdb** Files associated with the spatial databases.
* **viz** Directory containing ParaView Python scripts and other files for visualizing results.
* **output** Directory containing simulation output. It is created automatically when running the simulations.
* **README.md** README file containing a brief description of the various examples.

## Features Illustrated

{numref}`tab:example:reverse:2d:features` lists the features discussed in each of these 2-D reverse fault examples.

:::{admonition} TODO @brad
:class: error

Use script to create `features.md`.
:::

With the intent of illustrating features used in research simulations, we use HDF5 output and we make extensive use the most efficient implementations of spatial databases (UniformDB and ZeroDB).
We also use ParaView Python scripts for visualizing the output.
These scripts can be run within the ParaView GUI or outside the ParaView GUI, although the interaction is limited to rotating, translating, and zooming when run outside the ParaView GUI.


## Organization of Simulation Parameters

PyLith automatically reads in `pylithapp.cfg` from the current directory, if it exists.
As a result, we generally put all parameters common to a set of examples in this
file to avoid duplicating parameters across multiple files.
Because we often use a single mesh for multiple simulations in a directory, we place all
parameters related to our mesh and identifying the materials in our mesh in `pylithapp.cfg`.
We assign the bulk constitutive model and its parameters to each material in other files, because we generally vary those across the simulations.
In general, we place roller boundary conditions (Dirichlet boundary conditions constraining the degrees of freedom perpendicular to the boundary) on the lateral and bottom boundaries, so we include those in `pylithapp.cfg`.
In some simulations we will overwrite the values for parameters will values specific to a given example.
We also do the same thing for materials, since most of the examples use the default linear isotropic material.
This file is also a convenient place to put basic solver parameters and to turn on Pyre journals for displaying informational messages during a run/journalling debugging flags.

## Example Workflow

:::{toctree}
meshing.md
step01-gravity.md
step02-gravity-refstate.md
step03-gravity-incompressible.md
step04-surfload.md
step05-onefault.md
step06-twofaults-elastic.md
step07-twofaults-maxwell.md
step08-twofaults-powerlaw.md
:::
