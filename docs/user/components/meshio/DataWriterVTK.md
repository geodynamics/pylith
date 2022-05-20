# DataWriterVTK

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.DataWriterVTK`
:Journal name: `datawritervtk`

Writer of solution, auxiliary, and derived subfields to a VTK file.

Implements `DataWriter`.

## Pyre Properties

* `filename`=\<str\>: Name of VTK file.
  - **default value**: ''
  - **current value**: '', from {default}
* `float_precision`=\<int\>: Precision of floating point values in output.
  - **default value**: 6
  - **current value**: 6, from {default}
  - **validator**: (greater than 0)
* `time_constant`=\<dimensional\>: Values used to normalize time stamp in filename.
  - **default value**: 1*s
  - **current value**: 1*s, from {default}
  - **validator**: (greater than 0*s)
* `time_format`=\<str\>: C style format string for time stamp in filename.
  - **default value**: '%f'
  - **current value**: '%f', from {default}

## Example

Example of setting `DataWriterVTK` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[data_writer]
filename = domain_solution.vtk
time_format = %0.2f
time_constant = 1.0*year
float_precision = 6
:::

