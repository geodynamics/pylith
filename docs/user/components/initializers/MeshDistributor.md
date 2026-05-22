# MeshDistributor

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.initializers.MeshDistributor`
:Journal name: `mesh_refiner`

Distribute mesh.

Implements `InitializePhase`.

## Pyre Facilities

* `distributor`: Mesh distributor.
  - **current value**: 'mesh_distributor', from {default}
  - **configurable as**: mesh_distributor, distributor

## Example

Example of setting `MeshDistributor` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.mesh_initializer.phases.distribute_mesh]
distributor = pylith.topology.Distributor
:::

