# Initializer

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.initializers.Initializer`
:Journal name: `mesh_initializer`

Manager for reading and setting up a finite-element mesh.

## Pyre Facilities

* `phases`: Phases for mesh initialization.
  - **current value**: 'serial_phases', from {default}
  - **configurable as**: serial_phases, phases

## Example

Example of setting `Initializer` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.mesh_initializer]
phases = pylith.initializers.Serial

# Equivalent manual construction of Serial phases
phases = [read_mesh, reorder_mesh, distribute_mesh, insert_interfaces, refine_mesh]
read_mesh = pylith.initializers.MeshReader
reorder_mesh = pylith.initializers.MeshReordering
distribute_mesh = pylith.initializers.MeshDistributor
insert_interfaces = pylith.initializers.InsertInterfaces
refine_mesh = pylith.initializers.MeshRefiner
:::

