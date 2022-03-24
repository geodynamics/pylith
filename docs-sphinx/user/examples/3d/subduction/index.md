(sec:example:subduction:3d)=
# Examples for a 3D Subduction Zone

## Overview
This suite of examples demonstrates use of a wide variety of features and the general workflow often used in research simulations.
We base the model on the Cascadia subduction zone ({numref}`fig:example:subduction:3d:cascadia`).
These examples will focus on modeling the deformation associated with the the subducting slab, including interseismic deformation with aseismic slip (creep) and viscoelastic relaxation, coseismic slip on the slab interface and a splay fault, and slow slip events on the subduction interface.
We want to account for the 3-D material properties associated with different elastic properties for the subducting slab, mantle, continental crust, and an accretionary wedge.
To keep the computation time in these examples short, we limit our model to an 800 km {math}`\times` 800 km {math}`\times` 400 km domain and we will use a relatively coarse discretization.
For simplicity and to reduce complexity in constructing the mesh, we use a flat top surface (elevation of 0 with respect to mean sea level).

:::{figure-md} fig:example:subduction:3d:cascadia
<img src="figs/subduction3d_cascadia.*" alt="Cartoon of the Cascadia Subduction Zone showing the subduction of the Juan de Fuca Plate under the North American Plate. Source - [U.S. Geological Survey Fact Sheet 060-00](https://pubs.usgs.gov/fs/2000/fs060-00/)" width="100%"/>

Cartoon of the Cascadia Subduction Zone showing the subduction of the Juan de Fuca Plate under the North American Plate. Source - [U.S. Geological Survey Fact Sheet 060-00](https://pubs.usgs.gov/fs/2000/fs060-00/)
:::

{numref}`fig:example:subduction:3d:concept` shows our conceptual model with a slab, mantle, continental crust, and accretionary wedge.
We cut off the slab at a depth of 100 km.
We use a transverse geographic projection coordinate system with Portland, Oregon, as the origin in order to georeference our model.
In order to model the motion of the slab, we include a fault for the subduction interface (the interface between the top of the slab and the mantle, crust, and wedge), as well as a fault between the bottom of the slab and the mantle.

:::{figure-md} fig:example:subduction:3d:concept
<img src="figs/subduction3d_conceptualmodel.*" alt="Conceptual model based on the Cascadia Subduction Zone. The model includes the subduction slab (white), the mantle (green), continental crust (blue), and an accretionary wedge (red)." width="100%"/>

Conceptual model based on the Cascadia Subduction Zone. The model includes the subduction slab (white), the mantle (green), continental crust (blue), and an accretionary wedge (red).
:::

The files associated with this suite of examples are contained in the directory `examples/3d/subduction`.
This directory contains several subdirectories:

**mesh** Files used to construct the finite-element mesh using CUBIT/Trelis.

**spatialdb** Files associated with the spatial and temporal databases.

**viz** ParaView Python scripts and other files for visualizing results.

**output** Directory containing simulation output. It is created automatically when running the simulations.

## Features Illustrated

{numref}`tab:example:subduction:3d:features` lists the features discussed in each of these 3-D subduction zone examples.
 With the intent of illustrating features used in research simulations, we use HDF5 output and, we make extensive use the most efficient implementations of spatial databases (UniformDB and SimpleGridDB).
We also use ParaView Python scripts for visualizing the output.
These scripts can be run within the ParaView GUI or outside the ParaView GUI, although the interaction is limited to rotating, translating, and zooming when run outside the ParaView GUI.

:::{admonition} TODO
:class: error

Put in table: 3d_subduction_features.tex
:::

## Example Workflow

:::{toctree}
meshing.md
common-information.md
step01-axial-compression.md
step02-coseismic-slip.md
step03-aseismic-creep.md
step04-earthquake-cycle.md
step05-spontaneous-rupture.md
step06-prescribed-slow-slip.md
step07-inversion-greens-fns.md
step08-gravitational-stress-field.md
:::
