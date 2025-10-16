# MeshRefiner

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.initializers.MeshRefiner`
:Journal name: `mesh_refiner`

Refine mesh.

Implements `InitializePhase`.

## Pyre Facilities

* `refiner`: Mesh refiner.
  - **current value**: 'refineuniform', from {default}
  - **configurable as**: refineuniform, refiner

## Example

Example of setting `MeshRefiner` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.mesh_initializer.phases.refine_mesh]
refiner = pylith.topology.RefineUniform
:::

