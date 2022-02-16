(sec-user-run-pylith-utilities)=
# Utilities

The PyLith distribution contains several utilities for working with PyLith simulations and processing output.
These Python scripts are all installed into the same \filename{bin} directory as the `pylith` application with the exception of the \filename{pyre\_doc.py} script which is installed as part of Pythia/Pyre.

pyre_doc
: Display the Pyre properties and facilities available for a given component.

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

## pyre_doc.py

:::{note}
New in v3.0.0
:::

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
    normalizer=<component name>: Nondimensionalizer for problem.
        current value: 'nondimelasticquasistatic', from {default}
        configurable as: nondimelasticquasistatic, normalizer
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

## pylith_cfgsearch

:::{note}
New in v3.0.0
:::

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
caption: Example of running `pylith_cfgsearch` in `examples/2d/strikeslip`.
---
$ pylith_cfgsearch
step01_slip.cfg v1.0.0; requires PyLith >=3.0 and <4.0
    Coseismic prescribed slip with zero displacement Dirichlet boundary conditions.
    Authors: Brad Aagaard
    Keywords: example, 2D, strike slip, prescribed slip
    Features:
        Triangular cells, pylith.meshio.MeshIOCubit, pylith.problems.TimeDependent, pylith.materials.Elasticity,
        pylith.materials.IsotropicLinearElasticity, pylith.faults.FaultCohesiveKin, pylith.faults.KinSrcStep, field split
        preconditioner, Schur complement preconditioner, pylith.bc.DirichletTimeDependent, spatialdata.spatialdb.UniformDB,
        pylith.meshio.OutputSolnBoundary, pylith.meshio.DataWriterHDF5, Static simulation
    pylith step01_slip.cfg
step02_slip_velbc.cfg v1.0.0; requires PyLith >=3.0 and <4.0
    Coseismic prescribed slip with velocity Dirichlet boundary conditions.
    Authors: Brad Aagaard
    Keywords: example, 2D, strike slip, prescribed slip, velocity boundary conditions
    Features:
        Triangular cells, pylith.meshio.MeshIOCubit, pylith.problems.TimeDependent, pylith.materials.Elasticity,
        pylith.materials.IsotropicLinearElasticity, pylith.faults.FaultCohesiveKin, pylith.faults.KinSrcStep, field split
        preconditioner, Schur complement preconditioner, pylith.bc.DirichletTimeDependent, spatialdata.spatialdb.UniformDB,
        pylith.meshio.OutputSolnBoundary, pylith.meshio.DataWriterHDF5, Quasistatic simulation, spatialdata.spatialdb.SimpleDB
    pylith step02_slip_velbc.cfg
step03_multislip_velbc.cfg v1.0.0; requires PyLith >=3.0 and <4.0
    Coseismic prescribed slip with multiple ruptures and velocity Dirichlet boundary conditions.
    Authors: Brad Aagaard
    Keywords: example, 2D, strike slip, prescribed slip, multiple fault ruptures, velocity boundary conditions
    Features:
        Triangular cells, pylith.meshio.MeshIOCubit, pylith.problems.TimeDependent, pylith.materials.Elasticity,
        pylith.materials.IsotropicLinearElasticity, pylith.faults.FaultCohesiveKin, pylith.faults.KinSrcStep, field split
        preconditioner, Schur complement preconditioner, pylith.bc.DirichletTimeDependent, spatialdata.spatialdb.UniformDB,
        pylith.meshio.OutputSolnBoundary, pylith.meshio.DataWriterHDF5, Quasistatic simulation, spatialdata.spatialdb.SimpleDB
    pylith step03_multislip_velbc.cfg
```

```{code-block} console
---
caption: Example of running `pylith_cfgsearch` in `examples/2d/strikeslip`, limiting output to the description and keywords.
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
caption: Example of running `pylith_cfgsearch` in `examples/2d/strikeslip`, filtering search results to quasistatic simulations.
---
$ pylith_cfgsearch --features="Quasistatic simulation"
step02_slip_velbc.cfg v1.0.0; requires PyLith >=3.0 and <4.0
    Coseismic prescribed slip with velocity Dirichlet boundary conditions.
    Authors: Brad Aagaard
    Keywords: example, 2D, strike slip, prescribed slip, velocity boundary conditions
    Features:
        Triangular cells, pylith.meshio.MeshIOCubit, pylith.problems.TimeDependent, pylith.materials.Elasticity,
        pylith.materials.IsotropicLinearElasticity, pylith.faults.FaultCohesiveKin, pylith.faults.KinSrcStep, field split
        preconditioner, Schur complement preconditioner, pylith.bc.DirichletTimeDependent, spatialdata.spatialdb.UniformDB,
        pylith.meshio.OutputSolnBoundary, pylith.meshio.DataWriterHDF5, Quasistatic simulation, spatialdata.spatialdb.SimpleDB
    pylith step02_slip_velbc.cfg
step03_multislip_velbc.cfg v1.0.0; requires PyLith >=3.0 and <4.0
    Coseismic prescribed slip with multiple ruptures and velocity Dirichlet boundary conditions.
    Authors: Brad Aagaard
    Keywords: example, 2D, strike slip, prescribed slip, multiple fault ruptures, velocity boundary conditions
    Features:
        Triangular cells, pylith.meshio.MeshIOCubit, pylith.problems.TimeDependent, pylith.materials.Elasticity,
        pylith.materials.IsotropicLinearElasticity, pylith.faults.FaultCohesiveKin, pylith.faults.KinSrcStep, field split
        preconditioner, Schur complement preconditioner, pylith.bc.DirichletTimeDependent, spatialdata.spatialdb.UniformDB,
        pylith.meshio.OutputSolnBoundary, pylith.meshio.DataWriterHDF5, Quasistatic simulation, spatialdata.spatialdb.SimpleDB
    pylith step03_multislip_velbc.cfg
```

## pylith_runner

:::{note}
New in v3.0.0
:::

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
caption: Example of using `pylith_runner` to run all simulations in `examples/2d` (output omitted).
---
$ pylith_runner --path=examples/2d
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
The utility works with output from simulations with either prescribed slip and/or spontaneous rupture.
Currently, we compute the shear modulus from a user-specified spatial database at the centroid of the fault cells. In the future we plan to account for lateral variations in shear modulus across the fault when calculating the seismic moment.
The Python script is a Pyre application, so its parameters can be specified using `cfg` and command line arguments just like PyLith.

The Pyre properties and facilities include:

:output_filename: Filename for output of slip information.
:faults: Array of fault names.
:filename_pattern: Filename pattern in C/Python format for creating filename for each fault. Default is `output/fault_\%s.h5`.
:snapshots: Array of timestamps for slip snapshosts ([-1] means use last time step in file, which is the default).
:snapshot_units: Units for timestamps in array of snapshots.
:db_properties: Spatial database for elastic properties.
:coordsys: Coordinate system associated with mesh in simulation.

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

% End of file
