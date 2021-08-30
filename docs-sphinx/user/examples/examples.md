# Examples

## Overview

This chapter includes several suites of examples.
Each suite includes several "steps" which are examples that increase in complexity from one "step" to the next.
In some cases, a later step may make use of output from an earlier step; these cases are clearly documented.
{numref}`tab:examples:overview` classifies the level of difficulty of each example suite and provides a general description of the type of problems discussed.

:::{admonition} TODO
:class: error

Update obsolete section links.
:::

```{table} Overview of example suites.
:name: tab:examples:overview
| **Directory** |                                **Section(s)**                                 | **Difficulty** | **Description**                                                                                                            |
|:----------------|:-----------------------------------------------------------------------------:|:--------------:|:---------------------------------------------------------------------------------------------------------------------------|
| `twocells`      |   {ref}`sec:example:twotri3` - {ref}`sec:examples:twotet4-geoproj`    |     novice     | Toy problems with ASCII two-cell meshes.                                                                                         |
| `3d/hex8`        |                      {ref}`sec:example:3dhex8`                        |    beginner    | Illustration of most features using simple CUBIT box mesh.                                                                          |
| `3d/tet4`       |                      {ref}`sec:example:3dtet4`                        |    beginner    | Illustration of refinement using simple LaGriT box mesh.                                                                         |
| `bar_shearwave` | {ref}`sec:example:shearwave:tri3` - {ref}`sec:example:shearwave:hex8` |    beginner    | Illustration of wave propagation using simple shear beam.                                                                        |
| `2d/subduction` |                   {ref}`sec:example:subduction:2d`                    |  intermediate  | Illustration of coseismic, postseismic, and creep deformation using a 2-D subduction zone cross-section with a CUBIT mesh. |
| `2d/greensfns`  |                    {ref}`sec:example:greensfns2d`                     |  intermediate  | Illustration of computing static Green&rsquo;s functions for a strike-slip and reverse fault using a CUBIT mesh.                 |
| `3d/subduction` |                   {ref}`sec:example:subduction:3d`                    |  intermediate  | Illustration of most PyLith features for quasistatic deformation using a 3-D subduction zone with a CUBIT mesh.                      |
```

The `3d/subduction` example suite is the newest and most comprehensive.
Users wanting to use PyLith in their research should work through relevant beginner examples and then the `3d/subduction` examples.

### Prerequisites

Before you begin any of the examples, you will need to install PyLith following the instructions in {ref}`cha:installation`.
For more complex examples, you will also need either Trelis (available from <https://coreform.com/>), CUBIT (available to US federal government agencies from <https://cubit.sandia.gov/>) or LaGriT (available from <https://meshing.lanl.gov/>) mesh generation software to create the meshes.
If you do not wish to create your own mesh at this time, the meshes are also provided as part of the example.
The ParaView <https://www.paraview.org/> visualization package may be used to view simulation results.
ParaView 3 includes built-in documentation that is accessed by clicking on the Help menu item.
Some additional documentation is available on the ParaView Wiki site <https://www.paraview.org/Wiki/ParaView>.
You may use other visualization software, but some adaption from what is described here will be necessary.
Furthermore, you can complete a subset of the example using files provided (as described below), skipping the steps for which you do not have the proper software packages installed.

### Input Files

The files needed to work through the examples are found in the `examples` directory under the top-level PyLith directory.
There are five examples in `examples/twocells`, each consisting of just two cells (elements).
These very simple examples make use of PyLith mesh ASCII format to define the mesh.
This format is useful for understanding the basics of how PyLith works, since it is easy to create these files by hand.
More complex problems, such as those found in `examples/3d`, use external mesh generation software to create the meshes.
All of the files used in the example problems are extensively documented with comments.

(sec:ParaView:Python:scripts)=
## ParaView Python Scripts

**New in v2.2.1**

In some of the examples (currently only the 2D and 3D subduction zone examples) we provide ParaView Python scripts for visualizing the input finite-element mesh and the PyLith simulation results.
Some of these scripts are very generic and are easily reused; others are more specific to the examples.
The primary advantage of the ParaView Python scripts is that they make it easy to replicate visualizations, whether they are produced by the developers and regenerated by users.

There are several different ways to run the ParaView Python scripts:

* Within the ParaView GUI, select `Tools`&#8594;`Python Shell`.
Override the default parameters as desired (which we will discuss later in this section).
Click on the button, and navigate to the select the script you want to run.
* From a shell (terminal window) start ParaView from the command line with the `--script=FILENAME` where `FILENAME` is the relative or absolute path to the ParaView Python script.
Note that this method does not provide a mechanism for overriding the default parameters.
* Run the ParaView Python script directly from a shell (terminal window) via the command line.
You can use command line arguments to override the default values for the parameters.
If pvpython is not in your PATH, then you can run a script called `MY_SCRIPT` using: `PATH_TO_PVPYTHON/pvpython MY_SCRIPT.py`

:::{tip}
Running the ParaView Python script from within the ParaView GUI allows further manipulation of the data, which is not possible when running the ParaView Python script outside the ParaView GUI. When run outside the ParaView GUI, the interaction is limited to rotating, translating, and zooming.
:::

:::{important}
The ParaView Python scripts run Python via `pvpython`, which is a customized version of the Python interpreter included in the ParaView distribution. This is different from python provided with your operating system and/or the one included in the PyLith distribution. This means you cannot, in general, import Python modules provided with the PyLith distribution into ParaView.
:::

:::{tip}
In creating the ParaView Python scripts, we performed the steps within the GUI while capturing the commands using `Tools`&#8594;`Start Trace` and then `Tools`&#8594;`Start Trace`. This makes it very easy to create the Python script. Note that we have omitted superfluous commands in the trace when transferring the trace into a Python script. See the ParaView documentation for additional information about the Python API.
:::

### Overriding Default Parameters

We setup the ParaView Python scripts, so that when they are run from the command line in the main directory for a given example, e.g., `examples/3d/subduction`, the script will produce the output discussed in the manual.
If you start ParaView from the OS X Dock or a similar method, like a shortcut, then you will need to override at least the default values for the data file(s).

In order to override the default values from within the ParaView GUI, simply set the values within the Python shell.
For example, to set the value of the variable `EXODUS_FILE` to the absolute path of the input file,

```{code-block} python
---
caption: ParaView Python shell
---
>>> EXODUS_FILE = "/home/johndoe/pylith/examples/3d/subduction/mesh/mesh_tet.exo"
```

In this case, we use the Python os module to get the absolute path of the home directory and append the path to the Exodus file with the appropriate separators for the operating system.

:::{important}
In each of the ParaView Python scripts, the names of the variables and their default values are given by the DEFAULTS dictionary near the top of the file.
:::

## Examples

:::{toctree}
2d_box.md
3d_box.md
2d_strikeslip.md
2d_reverse.md
2d_subduction.md
3d_subduction.md
:::

## Additional Examples

### CUBIT Meshing Examples

The directory `examples/meshing` contains several examples of using CUBIT to construct finite-element meshes for complex geometry.
This includes features such as constructing nonplanar fault geometry from contours, constructing topography from a DEM, and merging sheet bodies (surfaces).
A separate examples discusses defining the discretization size using a vertex field in an Exodus-II file.
See the `README` files in the subdirectories for more detailed descriptions of these examples.

### Troubleshooting Examples

The directory `examples/troubleshooting` contains a few examples to practice troubleshooting a variety of user errors in parameters files and problem setup.
The files with the errors corrected are in `examples/troubleshooting/correct`.
Step-by-step corrections are discussed in the troubleshooting PyLith simulations sessions of the 2014, 2015, 2017, and 2019 PyLith tutorials (<https://wiki.geodynamics.org/software:pylith:start>).

### Code Verification Benchmarks

The CIG GitHub software repository <https://github.com/geodynamics/pylith_benchmarks> contains input files for a number of community benchmarks.
The benchmarks do not include the mesh files because they are so large; instead they include the CUBIT journal files that can be used to generate the meshes.
Most, but not all, of the input files in the repository are updated for PyLith v2.0.0, so you will need to modify them if you use another version of PyLith.
