# SolnDispPresTracStrainVelPdotTdotLagrange

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.SolnDispPresTracStrainVelPdotTdotLagrange`
:Journal name: `solndispprestracstrainvelpdottdotlagrange`

Python solution field with displacement, pressure, and trace strain subfields, along with their time derivatives, and a lagrange fault.

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
* `pressure_t`: PressureT subfield.
  - **current value**: 'subfieldpressure_t', from {default}
  - **configurable as**: subfieldpressure_t, pressure_t
* `trace_strain`: Trace strain subfield.
  - **current value**: 'subfieldtracestrain', from {default}
  - **configurable as**: subfieldtracestrain, trace_strain
* `trace_strain_t`: TraceStrainT subfield.
  - **current value**: 'subfieldtracestrain_t', from {default}
  - **configurable as**: subfieldtracestrain_t, trace_strain_t
* `velocity`: Velocity subfield.
  - **current value**: 'subfieldvelocity', from {default}
  - **configurable as**: subfieldvelocity, velocity

## Example

Example of setting `SolnDispPresTracStrainVelPdotTdotLagrange` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem]
solution = pylith.problems.SolnDispPresTracStrainVelPdotTdotLagrange
:::

