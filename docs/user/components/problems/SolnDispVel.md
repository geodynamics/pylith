# SolnDispVel

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.SolnDispVel`
:Journal name: `solndispvel`

Container for solution subfields with displacement and velocity subfields.

## Pyre Facilities

* `displacement`: Displacement subfield.
  - **current value**: 'subfielddisplacement', from {default}
  - **configurable as**: subfielddisplacement, displacement
* `velocity`: Velocity subfield.
  - **current value**: 'subfieldvelocity', from {default}
  - **configurable as**: subfieldvelocity, velocity

## Example

Example of setting `SolnDispVel` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem]
solution = pylith.problems.SolnDispVel
:::

