# Output

PyLith produces four kinds of output:

* Information written to stdout (default) or other devices that is controlled by the Pyre journal settings;
* Information written to stdout that is controlled by PETSc options;
* Information about the progress of the simulation written to a text file that is controlled by progress monitors (see {ref}`sec-user-progress-monitors`);
* Output of the solution over the domain, external boundary, or at discrete points (see {ref}`sec-user-solution-observers`); and
* Output of the solution or state variables over a material, interface, or external boundary (see {ref}`sec-user-physics-observers`).

(sec-user-progress-monitors)=
## Progress Monitors

The progress monitors make it easy to monitor the general progress of long simulations, especially on clusters where stdout is not always easily accessible.
The progress monitors update a simulation's current progress by writing information to a text file.
The information includes time stamps, percent completed, and an estimate of when the simulation will finish.

### `ProgressMonitorTime`

This is the default progress monitor for time-stepping problems.
The monitor calculates the percent completed based on the time at the current time step and the total simulated time of the simulation, not the total number of time steps (which may be unknown in simulations with adaptive time stepping).

:::{seealso}
[`ProgressMonitorTime` Component](../components/problems/ProgressMonitorTime.md)
:::

### `ProgressMonitorStep`

This is the default progress monitor for problems with a specified number of steps, such as Green's function problems.
The monitor calculates the percent completed based on the number of steps (e.g., Green's function impulses completed).

:::{seealso}
[`ProgressMonitorStep` Component](../components/problems/ProgressMonitorStep.md)
:::

## Observers

Observer objects manage transferring the solution and state variables to other objects, including output and external software.
Currently, the only type of observers implemented in PyLith are ones that produce output.

(sec-user-output-observers)=
### Output Observers

PyLith currently supports output to HDF5/Xdmf and VTK files, which can be imported directly into a number of visualization tools, such as ParaView, Visit, and MayaVi.
The HDF5 files can also be directly accessed via Matlab and Python modules, such as h5py.
PyLith v3.x supports output of the solution subfields, all auxiliary fields (material properties, boundary condition parameters, and fault interface parameters, etc.) and fields derived from the auxiliary field and/or solution, such as strain and stress.

The HDF5 writer provides parallel binary output, whereas the VTK writer provides serial ASCII output.
Additionally, with the VTK writer each time step is written to a separate file; the HDF5 writer puts all information for a given domain, boundary condition, or fault interface into a single file.

Output observers have a data writer (see {ref}`sec-user-data-writers`) for writing the data in a specified format and a trigger (see {ref}`sec-user-output-triggers`) for specifying how often to write the output in a time-dependent simulation.

:::{important}
Fields with a basis order of 1 are written as vertex fields, whereas fields with a basis order of 0 are written as cell fields.
Fields with a basis order of 0 are kept at a basis order of 0 when output.
Fields with a basis order of 1 or more can be output with a basis order of 0 or 1.
:::

*New in v4.0.0*

The output observers produce files with either diagnostic information (`info` files) or solution and state variable information (`data` files).
The default behavior is the files include all available information.
The `info` files include the auxiliary subfields at the beginning of the simulation along with  surface orientation information for faults and boundary conditions.
The `data` files include all solution subfields, state variables, and derived fields (fields computed from the solution, such as Cauchy stress and strain).

For boundary conditions the orientation information is provided in terms of x, y, and z components of the unit vectors for the surface normal and tangential directions.
In 3D the "vertical" tangential direction is the cross product of the surface normal and horizontal tangential direction; it is in the +z direction for a vertical boundary.
In the case of the horizontal boundary, the horizontal tangential direction is in the +x direction and the "vertical" tangential direction is in the +y direction.
For a fault surface the horizontal tangential direction generally corresponds to the along-strike direction and the "vertical" tangential direction generally corresponds to the up-dip direction; the exception is a 2D simulation for a vertical cross-section in which the "horizontal" tangential direction corresponds to the dip direction.

The orientation information is useful for transforming from components written in terms of a surface (for example, left lateral, reverse, and opening fault tractions), into the model coordinate system (for example, the global Cartesian coordinate system).
Given unit vector components $n_x$, $n_y$, $n_z$ (normal direction), $h_x$, $h_y$, $h_z$ (horizontal tangential direction or along-strike direction), and $v_x$, $v_y$, $v_z$ (vertical tangential direction or up-dip direction), we can transform a vector in the boundary (or fault) coordinate system ($a_n$, $a_h$, $a_v$) into the global coordinate system using

\begin{align}
a_x &= a_n n_x + a_h h_x + a_v v_x \\
a_y &= a_n n_y + a_h h_y + a_v v_y \\
a_z &= a_n n_z + a_h h_z + a_v v_z
\end{align}

This transformation is useful for plotting fault tractions and slip in 3D visualization tools such as ParaView that require vectors in the model coordinate system.

*New in v4.2.0*

Output observers can output information on a finer resolution mesh than the one used in the simulaiton.
This is useful when using a basis order of 2 or greater for the discretization of the solution and writing output at a basis order of 1.
In this case, output using a mesh that is 4-8 times finer and fields discretized using a basis order of 1 (usually required by visualization tools) will do a good job representing the fields with a basis order of 2 on the mesh used in the simulation.
This feature is used in [`examples/strikeslip-2d/step04-varslip`](../examples/strikeslip-2d/step04-varslip.md).

(sec-user-solution-observers)=
### Solution Observers

The solution observers get notified of updates to the solution.
{numref}`tab:solution:observers` lists the current implementations of solution observers, which are used for output.

```{table} Solution observers.
:name: tab:solution:observers
| Object               | Use Cases                                       |
| :------------------- | :---------------------------------------------- |
| `OutputSoln`         | Output of the solution over the domain          |
| `OutputSolnBoundary` | Output of the solutin over an external boundary |
| `OutputSolnPoints`   | Output of the solution at discrete points       |
```

:::{seealso}
[`OutputSoln` Component](../components/meshio/OutputSoln.md), [`OutputSolnBoundary` Component](../components/meshio/OutputSolnBoundary.md), and [`OutputSolnPoints` Component](../components/meshio/OutputSolnPoints.md).
:::

(sec-user-output-solution-points)=
#### Output at Discrete Points

In many situations with recorded observations, one would like to extract the solution at the same locations as the recorded observation.
=The locations are specified in a text file.

#### `PointsList` Reader

This object corresponds to a simple text file containing a list of points (one per line) where output is desired.
See {ref}`sec-user-file-formats-points-list` for file format specifications.
The points are specified in the coordinate system specified by `OutputSolnPoints`.
The coordinates will be transformed into the coordinate system of the mesh prior to interpolation. 

:::{seealso}
[`PointsList` Component](../components/meshio/PointsList.md)
:::

(sec-user-physics-observers)=
### Physics Observer

Analogous to the `OutputSoln` objects, which provide a means to output the solution, the physics objects (material, boundary conditions, and fault interfaces) have `OutputPhysics` objects to provide output of the solution, properties, state variables, etc.

(sec-user-data-writers)=
## Data Writers

(sec-user-data-writer-hdf5)=
### HDF5 Output

HDF5 files provide a flexible framework for storing simulation data with datasets in groups logically organized in a tree structure analogous to files in directories.
HDF5 output offers parallel, multi-dimensional array output in binary files, so it is much faster and more convenient than the VTK output which uses ASCII files and separate files for each time step.
Standards for organizing datasets and groups in HDF5 files do not exist for general finite-element software in geodynamics.
Consequently, PyLith uses its own simple layout show in {numref}`fig:hdf5:layout`.
In order for visualization tools, such as ParaView, to determine which datasets to read and where to find them in the hierarchy of groups within the HDF5 file, we create an Xdmf (eXtensible Data Model and Format, <https://www.xdmf.org> metadata file that provides this information.
This file is written when PyLith closes the HDF5 file at the end of the simulation.
In order to visualize the datasets in an HDF5 file, one simply opens the corresponding Xdmf file (the extension is `xmf`) in ParaView or Visit.
The Xdmf file contains the relative path to the HDF5 file so the files can be moved but must be located together in the same directory.

:::{important}
The Xdmf format supports representation of two- and three-dimensional coordinates of points, scalar fields, and three-dimensional vector and tensor fields but not two-dimensional vector or tensor fields.
Consequently, for two-dimensional vector fields we build a three-component vector from the two-component vector (x and y components) and a separate zero scalar field (z component).
For tensor fields, we create a scalar field for each of the tensor components, adding the component as a suffix to the name of the field.
:::

:::{figure-md} fig:hdf5:layout
<img src="figs/hdf5layout.*" alt="Diagram of HDF5 layout" width="100%"/>

General layout of a PyLith HDF5 file.
The orange rectangles with rounded corners identify the groups and the blue rectangles with sharp corners identify the datasets.
The dimensions of the data sets are shown in parentheses.
Most HDF5 files will contain either `vertex_fields` or `cell_fields` but not both.
:::

```{table} General ordering and names of vector and tensor components in HDF5 output.
:name: tab:output:components:order
| Vector Field Type| Components         |
|:-------------|:-----------------------|
| vector       | x, y, z                |
| tensor       | xx, yy, zz, xy, yz, xz |
```

PyLith provides two different data writers for HDF5 files.
The `DataWriterHDF5` object writes all information into the HDF5 file, whereas the `DataWriterHDF5Ext` object writes the data to external binary files and only the metadata to the HDF5 file.
HDF5 files do not contain self-correcting features that allow a file to be read if part of a dataset is corrupted.
This type of error can occur if a job terminates abnormally in the middle or at the end of a simulation on a large cluster or other parallel machine.
Fortunately, HDF5 also offers the ability to store datasets in external binary files with the locations specified by links in the HDF5 file.
Note that the use of external data files results in one data file per dataset in addition to the HDF5 and Xdmf files.
The external data files use the name of the HDF5 file with the dataset name added to the
prefix and the `h5` suffix replaced by `dat`.
The HDF5 files include relative paths to the external data files, so these files can also be moved, but they, too, must be kept together in the same directory.
This provides a more robust method of output because one can generate an HDF5 file associated with the uncorrupted portions of the external data files should an error occur.
Currently, PyLith does not include a utility to do this, but we plan to add one in a future release.
Thus, there are two options when writing PyLith output to HDF5 files: (1) including the datasets directly in the HDF5 files themselves using the `DataWriterHDF5` object or (2) storing the datasets in external binary files with just metadata in the HDF5 files using the `DataWriterHDF5Ext` object. Both methods provide similar performance because they will use MPI I/O if it is available.

:::{warning}
Storing the datasets within the HDF5 file in a parallel simulation requires that the HDF5 library be configured with the `--enable-parallel` option.
The binary PyLith packages include this feature and it is a default setting in building HDF5 via the PyLith Installer.
:::

Accessing the datasets for additional analysis or visualization is nearly identical in the two methods because the use of external data files is completely transparent to the user except for the presence of the additional files.
Note that in order for ParaView to find the HDF5 and external data files, it must be run from the same relative location where the simulation was run.
For example, if the simulation was run from a directory `work` and the HDF5/Xdmf files were written to `work/output`, then ParaView should be run from the `work` directory.

:::{seealso}
[`DataWriterHDF5` Component](../components/meshio/DataWriterHDF5.md) and [`DataWriterHDF5Ext` Component](../components/meshio/DataWriterHDF5Ext.md)
:::

#### HDF5 Utilities

HDF5 includes several utilities for examining the contents of HDF5 files.
`h5dump` is very handy for dumping the hierarchy, dimensions of datasets, attributes, and even the dataset values to stdout.

:::{code-block} console
# Dump the entire HDF5 file (not useful for large files).
$ h5dump mydata.h5

# Dump the hierarchy of an HDF5 file.
$ h5dump -n mydata.h5

# Dump the hierarchy with dataset dimensions and attributes.
$ h5dump -H mydata.h5

# Dump dataset 'vertices' in group '/geometry' to stdout.
$ h5dump -d /geometry/vertices mydata.h5
:::

We also provide a utility `pylith_genxdmf` (see {ref}`sec-user-pylith-genxdmf) that generates an appropriate Xdmf file from a PyLith HDF5 file.
This is very useful if you add fields to HDF5 files in post-processing and wish to view the results in ParaView or Visit.

(sec-user-data-writer-vtk)=
### VTK Output

PyLith writes VTU files for VTK output, which are in the XML format.
The XML files contain metadata in XML tags and raw data for the vertex coordinates, mesh topology, and fields over vertices and/or cells.
Each time step is written to a different file.
The time stamp is included in the filename with the decimal point removed.
This allows automatic generation of animations with many visualization packages that use VTK files.
The default time stamp is the time in seconds, but this can be changed using the normalization constant to give a time stamp in years, tens of years, or any other value.

:::{warning}
We discourage use of the VTK output as the files are not easily read from post-processing scripts.
We strongly recommend using HDF5 output instead, which is the default starting in PyLith v3.0.0.
:::

:::{seealso}
[`DataWriterVTK` Component](../components/meshio/DataWriterVTK.md)
:::

(sec-user-output-triggers)=
## Output Triggers

:::{note}
New in PyLith v3.0.0
:::

By default PyLith will write the requested output after every time step.
In many cases we prefer to save the solution, state variables, etc at a coarser temporal resolution.
`OutputTriggerStep` controls the decimation of the output by time step, and `OutputTriggerTime` controls the decimation of the output via elasped time.
For a constant time step these can be equivalent.

### Decimate by time step

`OutputTriggerStep` decimates the output by skipping a user-specified number of time steps.

:::{seealso}
[`OutputTriggerStep` Component](../components/meshio/OutputTriggerStep.md)
:::

### Decimate by time

`OutputTriggerTime` decimates the output by skipping a user-specified elasped time between output time slices.

:::{tip}
Due to roundoff error in determining the simulation time over many time steps, a simulation may occasionally skip writing output unexpectedly when using `OutputTriggerTime`.
The best workaround is to use an `elapsed_time` that is a fraction of the time step size smaller than the desired elapsed time, such as 0.9999*year instead of 1.0*year.
:::

:::{seealso}
[`OutputTriggerTime` Component](../components/meshio/OutputTriggerTime.md)
:::
