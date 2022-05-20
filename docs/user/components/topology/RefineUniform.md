# RefineUniform

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.topology.RefineUniform`
:Journal name: `refineuniform`

Uniform global mesh refinement in parallel.

Implements `MeshRefiner`.

## Pyre Properties

* `levels`=\<int\>: Number of refinement levels.
  - **default value**: 1
  - **current value**: 1, from {default}
  - **validator**: (greater than or equal to 1)

## Example

Example of setting `RefineUniform` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# Refine mesh twice to reduce size of cell edges by a factor of 4.
[pylithapp.mesh_generator.refiner]
levels = 2
:::

