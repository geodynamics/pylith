(sec-examples)=
# Examples

## Overview

This section includes several suites of examples available in the `examples` directory.
Each suite includes several "steps" which are examples that increase in complexity from one "step" to the next.
In some cases, a later step may make use of output from an earlier step; these cases are clearly documented.
{numref}`tab:examples:overview` classifies the level of difficulty of each example suite and provides a general description of the type of problems discussed.

```{table} Overview of example suites.
:name: tab:examples:overview

| **Example Suite**                                   | **Difficulty** | **Description**                                                                                                        |
| :-------------------------------------------------- | :------------: | :--------------------------------------------------------------------------------------------------------------------- |
| [`box-2d`](box-2d/index.md)                         |     novice     | Simple axial and shear deformation in static and quasi-static simulations in 2D box with a mesh in an ASCII text file. |
| [`box-3d`](box-3d/index.md)                         |     novice     | Same as `2d/box` but with a 3D box and a mesh from Gmsh or Cubit.                                                      |
| [`strikeslip-2d`](strikeslip-2d/index.md)           |    beginner    | Prescribed coseismic slip and multiple earthquake ruptures in 2D with a mesh from Gmsh or Cubit.                       |
| [`reverse-2d`](reverse-2d/index.md)                 |    beginner    | Gravity, surface loads, and prescribed coseismic slip on multiple reverse faults in 2D with a mesh from Gmsh or Cubit. |
| [`crustal-strikeslip-2d`](crustal-strikeslip-2d/index.md)           |  intermediate  | Coseismic deformation on 2D intersecting strike-slip faults with a mesh from Gmsh or Cubit. |
| [`crustal-strikeslip-3d`](crustal-strikeslip-3d/index.md)           |  intermediate  | Coseismic deformation on 3D intersecting strike-slip faults with a mesh from Gmsh or Cubit. |
| [`subduction-2d`](subduction-2d/index.md)           |  intermediate  | Coseismic, postseismic, and creep deformation using a 2D subduction zone cross-section with a mesh from Gmsh or Cubit. |
| [`subduction-3d`](subduction-3d/index.md)           |  intermediate  | Close to research-complexity for a 3D subduction zone with a mesh from Cubit.                                          |
| [`magma-2d`](magma-2d/index.md)                     |  intermediate  | Magma reservoir using poroelasticity.                                                                                  |
| [`poroelastic-outerrise-2d`](poroelastic-outerrise-2d/index.md)                 |  intermediate  | Hydration in the outer-rise of subducting oceanic lithosphere using poroelasticity.                                                                |
| [`troubleshooting-2d`](troubleshooting-2d/index.md) |     novice     | Troubleshooting errors in simulation in put files.                                                                     |
```

The `crustal-strikeslip-2d`, `crustal-strikeslip-3d`, and `subduction-3d` example suites are the most advanced.
If you want to use PyLith in your research, we recommend that you work through relevant beginner examples and then these more sophisticated examples.

:::{tip}
You can use the `pylith_cfgsearch` utility (see {ref}`sec-user-run-pylith-utilities`) to search for examples based on keywords and features.
:::

### Prerequisites

Before you begin any of the examples, you will need to install PyLith following the instructions in {ref}`sec-install`.
You should also read {ref}`sec-user-run-pylith`.
Complete sets of input files are included in the examples.
However, if you wish to generate the finite-element meshes yourself, you will also need Gmsh (available from <https://gmsh.info> and included in the PyLith binary package), Coreform Cubit (available from <https://coreform.com/>), or CUBIT (available to US federal government agencies from <https://cubit.sandia.gov/>).

### Input Files

The files needed to work through the examples are found in the `examples` directory under the top-level PyLith directory.
All of the files used in the example problems are extensively documented with comments.

### Visualizing PyLith Output

You can visualize the results using the {ref}`sec-user-run-pylith-viz` utility or other general scientific visualization tools, such as [ParaView](https://www.paraview.org/).
ParaView includes built-in documentation that is accessed by clicking on the Help menu item.
Some additional documentation is available on the ParaView Wiki site <https://www.paraview.org/Wiki/ParaView>.

## Examples

:::{toctree}
box-2d/index.md
box-3d/index.md
strikeslip-2d/index.md
reverse-2d/index.md
crustal-strikeslip-2d/index.md
crustal-strikeslip-3d/index.md
subduction-2d/index.md
subduction-3d/index.md
magma-2d/index.md
poroelastic-outerrise-2d/index.md
troubleshooting-2d/index.md
examples-other.md
:::
