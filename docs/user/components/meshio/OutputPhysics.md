# OutputPhysics

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.OutputPhysics`
:Journal name: `outputphysics`

Output for objects implementing physics (materials and boundary conditions).

:::{tip}
Most output information can be configured at the problem level using the [`ProblemDefaults` Component](../problems/ProblemDefaults.md).
:::

Implements `OutputObserver`.

## Pyre Facilities

* `trigger`: Trigger defining how often output is written.
  - **current value**: 'outputtriggerstep', from {default}
  - **configurable as**: outputtriggerstep, trigger
* `writer`: Writer for data.
  - **current value**: 'datawriterhdf5', from {default}
  - **configurable as**: datawriterhdf5, writer

## Pyre Properties

* `data_fields`=\<list\>: Names of solution, auxiliary, and derived subfields to include in data output.
  - **default value**: ['all']
  - **current value**: ['all'], from {default}
* `info_fields`=\<list\>: Names of auxiliary subfields to include in info output.
  - **default value**: ['all']
  - **current value**: ['all'], from {default}
* `output_basis_order`=\<int\>: Basis order for output.
  - **default value**: 1
  - **current value**: 1, from {default}
  - **validator**: (in [0, 1])

## Example

Example of setting `OutputPhysics` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[observer]
# Skip two time steps between output.
output_trigger = pylith.meshio.OutputTriggerStep
output_trigger.num_skip = 2

# Write output to HDF5 file with name `boundary_xpos.h5`.
writer = pylith.meshio.DataWriterHDF5
writer.filename = boundary_xpos.h5

output_basis_order = 1
:::

