# PyLith Application

The top-level object is the PyLith application with five facilities:

:metadata: Simulation metadata;
:petsc: PETSc settings.
:mesher: Importer for the finite-element mesh;
:problem: Problem to run, including materials, boundary conditions, etc.; and
:dump_parameters: Dump parameters used and version information to file.

The `mesher` facility handles getting the finite-element mesh information, such as importing a mesh from a file.
The `problem` facility defines of the boundary value problem to solve and contains the most information.
It defines the finite-element discretization, materials, faults, initial conditions, boundary conditions, and output.

:::{seealso}
[`PyLithApp` Component](../components/apps/PyLithApp.md)
:::

## Simulation Metadata

:::{note}
New in v3.0.0
:::

We use metadata to provide a concise summary of a simulation.
The metadata gives structure to information previously placed in comments within the parameter files while also making this information machine readable.
The metadata makes it possible for a Python script to launch an entire suite of simulations and search for simulation parameter files based on the metadata content (see {ref}`sec-user-run-pylith-utilities` for more information).
For example, users can search examples that use a given feature.

At a minumum the metadata must include: (1) a description, (2) the command line arguments necessary to run the simulation, and (3) the PyLith version(s) that are compatible with the input files.
We strongly encourage users to include all of the metadata in their own PyLith parameter files.

:::{seealso}
[`SimulationMetadata` Component](../components/utils/SimulationMetadata.md)
:::

% End of file
