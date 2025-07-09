# MeshImporter

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.topology.MeshImporter`
:Journal name: `meshimporter`

Base class for reading a finite-element mesh from files.

Implements `MeshGenerator`.

## Pyre Facilities

* `distributor`: Distributes mesh among processes.
  - **current value**: 'mesh_distributor', from {default}
  - **configurable as**: mesh_distributor, distributor
* `reader`: Reader for mesh file.
  - **current value**: 'meshioascii', from {default}
  - **configurable as**: meshioascii, reader
* `refiner`: Performs uniform global mesh refinement after distribution among processes (default is no refinement).
  - **current value**: 'refiner', from {default}
  - **configurable as**: refiner

## Pyre Properties

* `check_topology`=\<bool\>: Check topology of imported mesh.
  - **default value**: True
  - **current value**: True, from {default}
* `reorder_mesh`=\<bool\>: Reorder mesh using reverse Cuthill-McKee.
  - **default value**: True
  - **current value**: True, from {default}

## Example

Example of setting `MeshImporter` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.meshimporter]
reorder_mesh = True
check_topology = True
reader = pylith.meshio.MeshIOCubit
refiner = pylith.topology.RefineUniform
:::

