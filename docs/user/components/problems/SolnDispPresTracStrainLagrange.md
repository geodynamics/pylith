# SolnDispPresTracStrainLagrange

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.SolnDispPresTracStrainLagrange`
:Journal name: `SolnDispPresTracStrainLagrange`

Container for solution subfields with displacement, pore pressure, trace strain subfields, and fault Lagrange multiplier subfields.

## Pyre Facilities

* `displacement`: Displacement subfield.
  - **current value**: 'subfielddisplacement', from {default}
  - **configurable as**: subfielddisplacement, displacement
* `lagrange_multiplier_fault`: Fault Lagrange multiplier subfield.
  - **current value**: 'subfieldlagrangefault', from {default}
  - **configurable as**: subfieldlagrangefault, lagrange_multiplier_fault
* `pressure`: Pressure subfield.
  - **current value**: 'subfieldpressure', from {default}
  - **configurable as**: subfieldpressure, pressure
* `trace_strain`: Trace strain subfield.
  - **current value**: 'subfieldtracestrain', from {default}
  - **configurable as**: subfieldtracestrain, trace_strain

## Example

Example of setting `SolnDispPresTracStrainLagrange` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem]
solution = pylith.problems.SolnDispPresTracStrainLagrange
:::

