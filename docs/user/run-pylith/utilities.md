(sec-user-run-pylith-utilities)=
# Utilities

The PyLith distribution contains several utilities for working with PyLith simulations and processing output.
These Python scripts are all installed into the same `bin` directory as the `pylith` application with the exception of the `pyre_doc.py` script which is installed as part of Pythia/Pyre.

pylith_viz
: Visualize output from PyLith

pyre_doc
: Display the Pyre properties and facilities available for a given component.

pylith_convertmesh
: Convert an Exodus II file from Cubit or a Gmsh file to an HDF5 file in the PETSc mesh format. 

pylith_cfgsearch
: Search and display metadata in `.cfg` files.

pylith_runner
: Run all PyLith simulations in a given path.

pylith_dumpparameters
: Dump simulation parameters, including default values, to a file for viewing.

pylith_eqinfo
: Compute earthquake rupture metrics from PyLith output.

pylith_genxdmf
: Generate Xdmf files from HDF5 files written by PyLith.

pylith_powerlaw_gendb
: Generate a spatial database with power-law bulk rheology parameters for PyLith.

(sec-user-run-pylith-viz)=
## pylith_viz

*New in v4.1.0.*

This utility provides an interactive graphical user interface for plotting output from PyLith simulations using [PyVista](https://docs.pyvista.org/version/stable/).
We demonstrate use of this utility in the examples.

:::{important}
`pylith_viz` works with **HDF5 files** that use the same layout as those written by PyLith.
:::

```{code-block} bash
pylith_viz --filenames FILENAMES [--help] subcommand subcommand_options
```

:--help: Display help information for script.
:--filenames FILENAMES: Comma separated list of relative paths to HDF5 files from PyLith.
:subcommand: Subcommand to use for plotting data. Available subcommands are `plot_mesh`, `plot_field`, and `warp_grid`.
:subcommand_options: Options specific to subcommand.

### plot_mesh Subcommand

Plot the mesh and, optionally, the mesh quality.
When plotting the mesh quality, you can use the slider to adjust the threshold for the range of metric values shown.
This makes it easier to identify cells with the poorest and best mesh quality.

```{code-block} bash
pylith_viz --filenames FILENAMES plot_mesh [--mesh-quality METRIC] [--help]
```
:--help: Display help information for subcommand.
:--mesh-quality METRIC: Plot the mesh quality using METRIC. Refer to the [VTK mesh quality documentation](https://vtk.org/doc/nightly/html/classvtkMeshQuality.html) for definitions of the metrics.

### plot_field Subcommand

Plot a diagnostic, solution, or derived field.
If the data has more than one time step, then a slider will be enabled so that you can plot the data at different time steps.
Additionally, you can press the `n` key to increment time and the `p` key to decrement time.

The colormap and its range are selected based upon the component.
When you plot the `magnitude` component, the plotting uses the `plasma` sequential colormap with the minimum value set to 0 and the maximum value set to the maximum absolute magnitude over the time span of the data.
When you plot a vector or tensor component, the plotting uses the a blue-red diverging colormap with zero at the center of the colormap and the range set to the maximum absolute value of the component over the time span of the data.

```{code-block} bash
pylith_viz --filenames FILENAMES plot_field [--field FIELD_NAME] [--component COMPONENT] [--hide-edges] [--help]
```

:--help: Display help information for subcommand.
:--field FIELD_NAME: Plot field FIELD_NAME (default: `displacement`).
:--component COMPONENT: Component to display (default: `magnitude`).
:--hide-edges: Do not show cell edges.

### warp_grid Subcommand

Plot a solution field and warp the grid according to the displacement vector.
You can specify a scale factor for how much the deformation is exaggerated using a slider or a command line argument.
The maximum exaggeration is 5.0 times the initial exaggeration.

If the data has more than one time step, then a slider will be enabled so that you can plot the data at different time steps.
Additionally, you can press the `n` key to increment time and the `p` key to decrement time.

The colormap and its range are selected based upon the component.
When you plot the `magnitude` component, the plotting uses the `plasma` sequential colormap with the minimum value set to 0 and the maximum value set to the maximum absolute magnitude over the time span of the data.
When you plot a vector or tensor component, the plotting uses the a blue-red diverging colormap with zero at the center of the colormap and the range set to the maximum absolute value of the component over the time span of the data.

```{code-block} bash
pylith_viz --filenames FILENAMES warp_grid [--field FIELD_NAME] [--component COMPONENT] [--exaggeration EXAGGERATION] [--hide-outline] [--hide-edges] [--help]
```

:--help: Display help information for subcommand.
:--field FIELD_NAME: Plot field FIELD_NAME (default: `displacement`).
:--component COMPONENT: Component to display (default: `magnitude`).
:--exaggeration EXAGGERATION: Scale factor for exaggerating the displacement deformation (default: 1000).
:--hide-outline: Turn off showing the grid outline with the warped grid.
:--hide-edges: Do not show cell edges on warped grid.

## pyre_doc.py

*New in v3.0.0.*

This utility is part of the Pythia/Pyre framework.
It will be installed to the `bin` directory where Pythia/Pyre is installed.
This utility extracts the Python docstrings and help information for Pyre components.

:::{warning}
This utility does not work on the PyLithApp application object because it is a mpi.Application object.
:::

```{code-block} bash
pyre_doc.py [--help] [--short] OBJECT
```

:--help: Display help information for script.
:--short: Display only the docstrings for the specified module or object. The default is to display the information for the specified object and all child objects, such as classes within a module.
:OBJECT: Full path to Python module or object, such as `pylith.problems.Problem` (module) or `pylith.problems.Problem.Problem` (object).

```{code-block} console
---
caption: Example of running `pyre_doc.py` on `pylith.problems.Problem`
---
$ pyre_doc.py pylith.problems.Problem
No help available for module pylith.problems.Problem.

class Problem
Python abstract base class for crustal dynamics problems.

    FACTORY: problem

facilities of 'problem':
    bc=<component name>: Boundary conditions.
        current value: 'emptybin', from {default}
        configurable as: emptybin, bc
    defaults=<component name>: Default options for problem.
        current value: 'problem_defaults', from {default}
        configurable as: problem_defaults, defaults
    gravity_field=<component name>: Database used for gravity field.
        current value: 'nullcomponent', from {default}
        configurable as: nullcomponent, gravity_field
    interfaces=<component name>: Interior surfaces with constraints or constitutive models.
        current value: 'emptybin', from {default}
        configurable as: emptybin, interfaces
    materials=<component name>: Materials in problem.
        current value: 'homogeneous', from {default}
        configurable as: homogeneous, materials
    scales=<component name>: Nondimensionalizer for problem.
        current value: 'nondimelasticquasistatic', from {default}
        configurable as: nondimelasticquasistatic, scales
    solution=<component name>: Solution field for problem.
        current value: 'solution', from {default}
        configurable as: solution
    solution_observers=<component name>: Observers (e.g., output) for solution.
        current value: 'singlesolnobserver', from {default}
        configurable as: singlesolnobserver, solution_observers
properties of 'problem':
    formulation=<str>: Formulation for equations.
        default value: 'quasistatic'
        current value: 'quasistatic', from {default}
        validator: (in ['quasistatic', 'dynamic', 'dynamic_imex'])
    help=<bool>: prints a screen that describes my traits
        default value: False
        current value: False, from {default}
    help-components=<bool>: prints a screen that describes my subcomponents
        default value: False
        current value: False, from {default}
    help-persistence=<bool>: prints a screen that describes my persistent store
        default value: False
        current value: False, from {default}
    help-properties=<bool>: prints a screen that describes my properties
        default value: False
        current value: False, from {default}
    solver=<str>: Type of solver to use ['linear', 'nonlinear'].
        default value: 'linear'
        current value: 'linear', from {default}
        validator: (in ['linear', 'nonlinear'])
```

(sec-user-run-pylith-convertmesh)=
## pylith_convertmesh

*New in v5.0.0.*

This utility converts a finite-element mesh from one file format to another.
The primary use case is converting an Exodus II file from Cubit or a Gmsh file to an HDF5 in the PETSc mesh format.
This is a Pyre application, so you can get help using the same command line tools as those for PyLith.
The default reader and writer are `MeshIOPetsc` components, so you do not need to change the reader or writer component to convert a Gmsh file to an HDF5 in the PETSc mesh format.

```{code-block} console
---
caption: Example of running `pylith_convertmesh` with command line arguments for getting help.
---
$ pylith_convertmesh [--help] [--help-components] [--help-properties]
```

```{code-block} console
---
caption: Example of running `pylith_convertmesh` to convert a Gmsh file `mesh_tri.msh` to an HDF5 file in PETSc mesh format `mesh_tri.h5`.
---
$ pylith_convertmesh --reader.filename=mesh_tri.mesh --writer.filename=mesh_tri.h5
```

```{code-block} console
---
caption: Example of running `pylith_convertmesh` to convert an Exodus II file `mesh_tri.exo` to an HDF5 file in PETSc mesh format `mesh_tri.h5`.
---
$ pylith_convertmesh --reader=pylith.meshio.MeshIOCubit --reader.filename=mesh_tri.exo --writer.filename=mesh_tri.h5
```

## pylith_cfgsearch

*New in v3.0.0.*

This utility searches and displays the metadata in `.cfg` files based on criteria provided via the command line.

```{code-block} bash
pylith_cfgsearch [--help] [--path SEARCHPATH] [--display DISPLAY] [--verbose] [--keywords KEYWORDS]
    [--features FEATURES] [--authors AUTHORS] [--version VERSION] [--incompatible]
```

:--help: Display help information for script.
:--path SEARCHPATH: Search path for `.cfg` files (default: `./`).
:--display DISPLAY: List of metadata to display in search results (default: all).
:--keywords KEYWORDS: Comma delimited list of keywords for filtering search results (default: None).
:--features FEATURES: Comma delimited list of features for filtering search results (default: None).
:--authors AUTHORS: Comma delimited list of authors for filtering search results (default: None).
:--version VERSION: PyLith version for filtering search results (default: None).
:--verbose: Report missing metadata (default: False).
:--incompatible: Filter search results to show incompatible parameter files (default: False).

```{code-block} console
---
caption: Example of running `pylith_cfgsearch` in `examples/strikeslip-2d`.
---
$ pylith_cfgsearch
step01_slip.cfg v1.0.0; requires PyLith >=3.0 and <5.0
    Coseismic prescribed slip with zero displacement Dirichlet boundary conditions.
    Authors: Brad Aagaard
    Keywords: example, 2D, strike slip, prescribed slip
    Features:
        Triangular cells, pylith.meshio.MeshIOCubit, pylith.problems.TimeDependent, pylith.materials.Elasticity,
        pylith.materials.IsotropicLinearElasticity, pylith.faults.FaultCohesiveKin, pylith.faults.KinSrcStep, field split
        preconditioner, Schur complement preconditioner, pylith.bc.DirichletTimeDependent, spatialdata.spatialdb.UniformDB,
        pylith.meshio.OutputSolnBoundary, pylith.meshio.DataWriterHDF5, Static simulation
    pylith step01_slip.cfg
step02_slip_velbc.cfg v1.0.0; requires PyLith >=3.0 and <5.0
    Coseismic prescribed slip with velocity Dirichlet boundary conditions.
    Authors: Brad Aagaard
    Keywords: example, 2D, strike slip, prescribed slip, velocity boundary conditions
    Features:
        Triangular cells, pylith.meshio.MeshIOCubit, pylith.problems.TimeDependent, pylith.materials.Elasticity,
        pylith.materials.IsotropicLinearElasticity, pylith.faults.FaultCohesiveKin, pylith.faults.KinSrcStep, field split
        preconditioner, Schur complement preconditioner, pylith.bc.DirichletTimeDependent, spatialdata.spatialdb.UniformDB,
        pylith.meshio.OutputSolnBoundary, pylith.meshio.DataWriterHDF5, Quasi-static simulation, spatialdata.spatialdb.SimpleDB
    pylith step02_slip_velbc.cfg
step03_multislip_velbc.cfg v1.0.0; requires PyLith >=3.0 and <5.0
    Coseismic prescribed slip with multiple ruptures and velocity Dirichlet boundary conditions.
    Authors: Brad Aagaard
    Keywords: example, 2D, strike slip, prescribed slip, multiple fault ruptures, velocity boundary conditions
    Features:
        Triangular cells, pylith.meshio.MeshIOCubit, pylith.problems.TimeDependent, pylith.materials.Elasticity,
        pylith.materials.IsotropicLinearElasticity, pylith.faults.FaultCohesiveKin, pylith.faults.KinSrcStep, field split
        preconditioner, Schur complement preconditioner, pylith.bc.DirichletTimeDependent, spatialdata.spatialdb.UniformDB,
        pylith.meshio.OutputSolnBoundary, pylith.meshio.DataWriterHDF5, Quasi-static simulation, spatialdata.spatialdb.SimpleDB
    pylith step03_multislip_velbc.cfg
```

```{code-block} console
---
caption: Example of running `pylith_cfgsearch` in `examples/strikeslip-2d`, limiting output to the description and keywords.
---
$ pylith_cfgsearch --display=description,keywords
step01_slip.cfg
    Coseismic prescribed slip with zero displacement Dirichlet boundary conditions.
    Keywords: example, 2D, strike slip, prescribed slip
step02_slip_velbc.cfg
    Coseismic prescribed slip with velocity Dirichlet boundary conditions.
    Keywords: example, 2D, strike slip, prescribed slip, velocity boundary conditions
step03_multislip_velbc.cfg
    Coseismic prescribed slip with multiple ruptures and velocity Dirichlet boundary conditions.
    Keywords: example, 2D, strike slip, prescribed slip, multiple fault ruptures, velocity boundary conditions
```

```{code-block} console
---
caption: Example of running `pylith_cfgsearch` in `examples/strikeslip-2d`, filtering search results to quasistatic simulations.
---
$ pylith_cfgsearch --features="Quasi-static simulation"
step02_slip_velbc.cfg v1.0.0; requires PyLith >=3.0 and <5.0
    Coseismic prescribed slip with velocity Dirichlet boundary conditions.
    Authors: Brad Aagaard
    Keywords: example, 2D, strike slip, prescribed slip, velocity boundary conditions
    Features:
        Triangular cells, pylith.meshio.MeshIOCubit, pylith.problems.TimeDependent, pylith.materials.Elasticity,
        pylith.materials.IsotropicLinearElasticity, pylith.faults.FaultCohesiveKin, pylith.faults.KinSrcStep, field split
        preconditioner, Schur complement preconditioner, pylith.bc.DirichletTimeDependent, spatialdata.spatialdb.UniformDB,
        pylith.meshio.OutputSolnBoundary, pylith.meshio.DataWriterHDF5, Quasi-static simulation, spatialdata.spatialdb.SimpleDB
    pylith step02_slip_velbc.cfg
step03_multislip_velbc.cfg v1.0.0; requires PyLith >=3.0 and <5.0
    Coseismic prescribed slip with multiple ruptures and velocity Dirichlet boundary conditions.
    Authors: Brad Aagaard
    Keywords: example, 2D, strike slip, prescribed slip, multiple fault ruptures, velocity boundary conditions
    Features:
        Triangular cells, pylith.meshio.MeshIOCubit, pylith.problems.TimeDependent, pylith.materials.Elasticity,
        pylith.materials.IsotropicLinearElasticity, pylith.faults.FaultCohesiveKin, pylith.faults.KinSrcStep, field split
        preconditioner, Schur complement preconditioner, pylith.bc.DirichletTimeDependent, spatialdata.spatialdb.UniformDB,
        pylith.meshio.OutputSolnBoundary, pylith.meshio.DataWriterHDF5, Quasi-static simulation, spatialdata.spatialdb.SimpleDB
    pylith step03_multislip_velbc.cfg
```

## pylith_runner

*New in v3.0.0.*

The runner utility searches a directory path for `.cfg` files with `arguments` in the simulation metadata (see {ref}`sec-user-run-pylith-pylithapp` for details).
It uses those arguments to run PyLith simulations.

:::{important}
The runner utility only knows how to run PyLith simulations, it does not know how to do pre- or post-processing.
:::

```{code-block} bash
pylith_runner [--help] [--path SEARCHPATH] [--verbose]
```

:--help: Display help information for script.
:--path SEARCHPATH: Search path for `.cfg` files (default: `./`).
:--verbose: Report missing metadata (default: False).

```{code-block} console
---
caption: Example of using `pylith_runner` to run all simulations in `examples/box-2d` (output omitted).
---
$ pylith_runner --path=examples/box-2d
```

(sec-user-run-pylith-pylith-dumpparameters)=
## pylith_dumpparameters

This utility dumps all PyLith parameters collected from `cfg` files, the command line, and any other sources, into a file.
The output file can be a simple text file or a JSON file.
The JSON file can be viewed using the [PyLith Parameter Viewer](sec-user-run-pylith-parameter-gui).

```{code-block} bash
pylith_dumpparameters [--quiet] [--format=ascii|json] [--filnemame=FILENAME]
```

:--quiet: Don't include description, location, or aliases in output.
:--format=ascii|json: Format of output file (default=`json`).
:--filename=FILENAME: Name of output file.

## pylith_eqinfo

This utility computes the moment magnitude, seismic moment, seismic potency, and average slip at user-specified time snapshots from PyLith fault HDF5 output.
The utility works with output from simulations with either prescribed slip or spontaneous rupture.
The moment magnitudes, seismic moment, and seismic potentency for 2D simulation have limited value, because they use on a 1D fault.
Currently, we compute the shear modulus from a user-specified spatial database at the centroid of the fault cells.
In the future we plan to automatically account for lateral variations in shear modulus across the fault when calculating the seismic moment.
The average slip, rupture area, seismic moment, and seismic potency are all given in SI units.

```{tip}
From {cite:t}`Wu:Chen:2003` the seismic moment, $M_0$, when considering a laterial variation in the shear modulus across the fault is

\begin{equation}
M_0 = \mu_\mathit{eff} A D
\end{equation}

where $A$ is the rupture area, $D$ is the average slip, and $\mu_\mathit{eff}$ is the effective shear modulus given by

\begin{equation}
\mu_\mathit{eff} = \frac{1}{2} \left( \mu^+ + \mu^- \right) \left(1 - \left(\frac{\mu^+ - \mu^-}{\mu^+ + \mu^-} \right)^2 \right)
\end{equation}

where $\mu^+$ and $\mu^-$ are the shear modulus on each side of the fault.
```

The Python script is a Pyre application, so its parameters can be specified using `cfg` files (`eqinfoapp.cfg` will be read if it exists) or command line arguments just like PyLith.

The Pyre properties and facilities include:

:output_filename: Filename for output of slip information.
:coordsys: Coordinate system associated with mesh in simulation.
:faults: Array of fault names.
:filename_pattern: Filename pattern in C/Python format for creating filename for each fault. Default is `output/fault_\%s.h5`.
:snapshots: Array of timestamps for slip snapshosts ([-1] means use last time step in file, which is the default).
:snapshot_units: Units for timestamps in array of snapshots.
:db_properties: Spatial database for elastic properties.
```{code-block} cfg
---
caption: Example `cfg` file for running `pylith_eqinfo` in `examples/strikeslip-2d` to compute the earthquake magnitude, seismic moment, and average slip for step04_varslip.
---
[eqinfoapp]
output_filename = output/step04_varslip-eqinfo.py

coordsys.space_dim = 2

faults = [fault]
filename_pattern = output/step04_varslip-%s.h5

db_properties = spatialdata.spatialdb.UniformDB
db_properties.description = Fault properties
db_properties.values = [density, Vs]
db_properties.data = [2500.0*kg/m**3, 3.46*km/s]
```

(sec-user-run-pylith-pylith-genxdmf)=
## pylith_genxdmf

This utility generates Xdmf files from HDF5 files that conform to the layout used by PyLith.
It is a simple Python script with a single command line argument with the file pattern of HDF5 files for which Xdmf files should be generated.
Typically, it is used to regenerate Xdmf files that get corrupted or lost due to renaming and moving.
It is also useful in updating Xdmf files when users add fields to HDF5 files during post-processing.

```{code-block} bash
pylith_genxdmf --files=FILE_OR_FILE_PATTERN
```

The default value for `FILE_OR_FILE_PATTERN` is `*.h5`.

:::{warning}
If the HDF5 files contain external datasets, then this utility should be run from the same relative path to the HDF5 files as when they were created.
For example, if a PyLith simulation was run from directory `work` and HDF5 files were generated in `output/work`, then the utility should be run from the directory `work`.
Furthermore, a visualization tool, such as ParaView, should also be started from the working directory `work`.
:::

(sec-user-run-pylith-pylith-powerlaw-gendb)=
## pylith_powerlaw_gendb

This Pyre application generates a `SimpleDB` spatial database file with power-law viscoelastic material properties for use with PyLith.
The inputs are spatial databases with values often available from laboratory experiments, such as activation energy, temperature, power-law coefficient, and the power-law exponent.
An additional parameter defines the units of the activation energy.
You must also specify either a reference stress or a reference strain rate.

You place all of the application parameters in `powerlaw_gendb.cfg`, which the application will read by default.
See {ref}`sec-user-examples-reverse-2d-step08` for an example of how to use `pylith_powerlaw_gendb`.

% End of file
