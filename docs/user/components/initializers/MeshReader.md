# MeshReader

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.initializers.MeshReader`
:Journal name: `mesh_reader`

Read mesh.

Implements `InitializePhase`.

## Pyre Facilities

* `reader`: Mesh reader.
  - **current value**: 'meshiopetsc', from {default}
  - **configurable as**: meshiopetsc, reader

## Example

Example of setting `MeshReader` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.mesh_initializer.phases.read_mesh]
reader = pylith.meshio.MeshIOPetsc
:::

