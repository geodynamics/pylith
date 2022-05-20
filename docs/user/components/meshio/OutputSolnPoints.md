# OutputSolnPoints

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.OutputSolnPoints`
:Journal name: `outputsolnpoints`

Output of solution subfields at discrete points in the domain.

:::{tip}
Most output information can be configured at the problem level using the [`ProblemDefaults` Component](../problems/ProblemDefaults.md).
:::

Implements `OutputSoln`.

## Pyre Facilities

* `reader`: Reader for points list.
  - **current value**: 'pointslist', from {default}
  - **configurable as**: pointslist, reader
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
* `label`=\<str\>: Label identifier for points (used in constructing default filenames).
  - **default value**: 'points'
  - **current value**: 'points', from {default}
* `output_basis_order`=\<int\>: Basis order for output.
  - **default value**: 1
  - **current value**: 1, from {default}
  - **validator**: (in [0, 1])

## Example

Example of setting `OutputSolnPoints` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[observer]
label = stations
data_fields = [displacement]

# List of points where we want output.
reader = pylith.meshio.PointsList
reader.filename = stations.txt

# Skip two time steps between output.
output_trigger = pylith.meshio.OutputTriggerStep
output_trigger.num_skip = 2

# Write output to HDF5 file with name `domain.h5`.
writer = pylith.meshio.DataWriterHDF5
writer.filename = domain.h5

output_basis_order = 1
:::

