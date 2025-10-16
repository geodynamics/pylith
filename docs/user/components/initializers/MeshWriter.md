# MeshWriter

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.initializers.MeshWriter`
:Journal name: `mesh_writer`

Write mesh.

Implements `InitializePhase`.

## Pyre Facilities

* `writer`: Mesh writer.
  - **current value**: 'meshiopetsc', from {default}
  - **configurable as**: meshiopetsc, writer

## Example

Example of setting `MeshWriter` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.mesh_initializer.phases.write_mesh]
writer = pylith.meshio.MeshIOPetsc
:::

