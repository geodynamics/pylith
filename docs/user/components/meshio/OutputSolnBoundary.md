# OutputSolnBoundary

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.OutputSolnBoundary`
:Journal name: `outputsolnsubset`

Output of solution subfields over an external boundary.

:::{tip}
Most output information can be configured at the problem level using the [`ProblemDefaults` Component](../problems/ProblemDefaults.md).
:::

Implements `OutputSoln`.

## Pyre Facilities

* `trigger`: Trigger defining how often output is written.
  - **current value**: 'outputtriggerstep', from {default}
  - **configurable as**: outputtriggerstep, trigger
* `writer`: Writer for data.
  - **current value**: 'datawriterhdf5', from {default}
  - **configurable as**: datawriterhdf5, writer

## Pyre Properties

* `data_fields`=\<list\>: Names of solution subfields to include in output.
  - **default value**: ['all']
  - **current value**: ['all'], from {default}
* `label`=\<str\>: Name of label identifier for external boundary.
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateLabel at 0x11f366af0>
* `label_value`=\<int\>: Value of label identifier for external boundary (tag of physical group in Gmsh files).
  - **default value**: 1
  - **current value**: 1, from {default}
* `output_basis_order`=\<int\>: Basis order for output.
  - **default value**: 1
  - **current value**: 1, from {default}
  - **validator**: (in [0, 1])

## Example

Example of setting `OutputSolnBoundary` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[observer]
data_fields = [displacement]

label = boundary_xpos

# Skip two time steps between output.
output_trigger = pylith.meshio.OutputTriggerStep
output_trigger.num_skip = 2

# Write output to HDF5 file with name `boundary_xpos.h5`.
writer = pylith.meshio.DataWriterHDF5
writer.filename = boundary_xpos.h5

output_basis_order = 1
:::

