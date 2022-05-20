# SolnDispPresTracStrainVelPdotTdot

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.SolnDispPresTracStrainVelPdotTdot`
:Journal name: `solndispprestracstrainveltdotpdot`

Container for solution subfields with displacement, pore pressure, and trace strain subfields, along with their time derivatices.

## Pyre Facilities

* `displacement`: Displacement subfield.
  - **current value**: 'subfielddisplacement', from {default}
  - **configurable as**: subfielddisplacement, displacement
* `pressure`: Pressure subfield.
  - **current value**: 'subfieldpressure', from {default}
  - **configurable as**: subfieldpressure, pressure
* `pressure_t`: Pressure_t subfield.
  - **current value**: 'subfieldpressure_t', from {default}
  - **configurable as**: subfieldpressure_t, pressure_t
* `trace_strain`: Trace strain subfield.
  - **current value**: 'subfieldtracestrain', from {default}
  - **configurable as**: subfieldtracestrain, trace_strain
* `trace_strain_t`: Trace strain_t subfield.
  - **current value**: 'subfieldtracestrain_t', from {default}
  - **configurable as**: subfieldtracestrain_t, trace_strain_t
* `velocity`: Velocity subfield.
  - **current value**: 'subfieldvelocity', from {default}
  - **configurable as**: subfieldvelocity, velocity

## Example

Example of setting `SolnDispPresTracStrainVelPdotTdot` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem]
solution = pylith.problems.SolnDispPresTracStrainVelPdotTdot
:::

