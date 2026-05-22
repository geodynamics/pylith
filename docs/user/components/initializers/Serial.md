# Serial

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.initializers.Serial`
:Journal name: `serial_phases`

Mesh initialization phases for reading in serial.

:::{seealso}
See [`Initializer` Component](Initializer.md).
:::

## Pyre Facilities

* `distribute_mesh`: Distribute mesh among processes.
  - **current value**: 'mesh_refiner', from {default}
  - **configurable as**: mesh_refiner, distribute_mesh
* `insert_interfaces`: Insert interfaces using PETSc mesh transform operation.
  - **current value**: 'mesh_insert_interfaces', from {default}
  - **configurable as**: mesh_insert_interfaces, insert_interfaces
* `read_mesh`: Read mesh in serial.
  - **current value**: 'mesh_reader', from {default}
  - **configurable as**: mesh_reader, read_mesh
* `refine_mesh`: Refine mesh.
  - **current value**: 'mesh_refiner', from {default}
  - **configurable as**: mesh_refiner, refine_mesh
* `reorder_mesh`: Reorder mesh using reverse Cuthill-McKee algorithm.
  - **current value**: 'reorder_mesh', from {default}
  - **configurable as**: reorder_mesh

## Example

Example of setting `Serial` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# Equivalent manual construction
phases = [read_mesh, distribute_mesh, insert_interfaces, refine_mesh]
read_mesh = pylith.initializers.MeshReader
distribute_mesh = pylith.initializers.MeshDistributor
insert_interfaces = pylith.initializers.InsertInterfaces
refine_mesh = pylith.initializers.MeshRefiner
:::

