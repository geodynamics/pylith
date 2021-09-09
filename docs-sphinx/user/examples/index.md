(cha:examples)=
# Examples

## Overview

This section includes several suites of examples.
Each suite includes several "steps" which are examples that increase in complexity from one "step" to the next.
In some cases, a later step may make use of output from an earlier step; these cases are clearly documented.
{numref}`tab:examples:overview` classifies the level of difficulty of each example suite and provides a general description of the type of problems discussed.

```{table} Overview of example suites.
:name: tab:examples:overview

| **Example Suite**                         | **Difficulty** | **Description** |
|:------------------------------------------|:------------:|:---------------------------------------------------------------------------------------------|
| [`2d/box`](2d/box/index.md)               | novice       | Simple axial and shear deformation in static and quasistatic simulations in 2D box with a CUBIT mesh.
| [`3d/box`](3d/box/index.md)               | novice       | Same as `2d/box` but with a 3D box with a CUBIT mesh.
| [`2d/strikeslip`](2d/strikeslip/index.md) | beginner     | Prescribed coseismic slip and multiple earthquake ruptures in 2D with a CUBIT mesh.
| [`2d/reverse`](2d/reverse/index.md)       | beginner     | Gravity, surface loads, and prescribed coseismic slip on multiple reverse faults in 2D with a CUBIT mesh.
| [`2d/subduction`](2d/subduction/index.md) | intermediate | Coseismic, postseismic, and creep deformation using a 2D subduction zone cross-section with a CUBIT mesh.
| [`3d/subduction`](3d/subduction/index.md) | intermediate | Close to research-complexity for a 3D subduction zone with a CUBIT mesh.
```

The `3d/subduction` example suite is the newest and most comprehensive.
Users wanting to use PyLith in their research should work through relevant beginner examples and then the `3d/subduction` examples.

### Prerequisites

Before you begin any of the examples, you will need to install PyLith following the instructions in {ref}`cha:installation`.
For more complex examples, you will also need either Coreform Cubit (available from <https://coreform.com/>), CUBIT (available to US federal government agencies from <https://cubit.sandia.gov/>) or LaGriT (available from <https://meshing.lanl.gov/>) mesh generation software to create the meshes.
If you do not wish to create your own mesh at this time, the meshes are also provided as part of the example.
The ParaView <https://www.paraview.org/> visualization package may be used to view simulation results.
ParaView 3 includes built-in documentation that is accessed by clicking on the Help menu item.
Some additional documentation is available on the ParaView Wiki site <https://www.paraview.org/Wiki/ParaView>.
You may use other visualization software, but some adaption from what is described here will be necessary.
Furthermore, you can complete a subset of the example using files provided (as described below), skipping the steps for which you do not have the proper software packages installed.

### Input Files

The files needed to work through the examples are found in the `examples` directory under the top-level PyLith directory.
All of the files used in the example problems are extensively documented with comments.

### Visualizing PyLith Output

See [ParaView Python Scripts](paraview-python.md) for a description of how to make use of the provided Python scripts for visualizing simulation output with ParaView.
Alternatively, you can use manually construct the visuzliation pipeline in several open-source visualization tools, such as ParaView and Visit. 

## Examples

:::{toctree}
2d/box/index.md
3d/box/index.md
2d/strikeslip/index.md
2d/reverse/index.md
2d/subduction/index.md
3d/subduction/index.md
examples-other.md
:::
