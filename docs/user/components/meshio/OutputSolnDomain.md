# OutputSolnDomain

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.OutputSolnDomain`
:Journal name: `outputsolndomain`

Output of solution subfields over the simulation domain.

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
* `output_basis_order`=\<int\>: Basis order for output.
  - **default value**: 1
  - **current value**: 1, from {default}
  - **validator**: (in [0, 1])
* `refine_levels`=\<int\>: Number of mesh refinement levels for output.
  - **default value**: 0
  - **current value**: 0, from {default}
  - **validator**: (greater than or equal to 0)

## Example

Example of setting `OutputSolnDomain` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[observer]
data_fields = [displacement]

# Skip two time steps between output.
trigger = pylith.meshio.OutputTriggerStep
trigger.num_skip = 2

# Write output to HDF5 file with name `domain.h5`.
writer = pylith.meshio.DataWriterHDF5
writer.filename = domain.h5

# Output with a basis order of 1 and refine mesh 3x (cells are 1/8 size).
# Refining the output mesh is useful with a basis order of 2 or greater in the solution.
output_basis_order = 1
refine_levels = 3
:::

