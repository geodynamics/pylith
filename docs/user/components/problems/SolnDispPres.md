# SolnDispPres

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.SolnDispPres`
:Journal name: `solndisppres`

Container for solution subfields with displacement and pressure subfields.

## Pyre Facilities

* `displacement`: Displacement subfield.
  - **current value**: 'subfielddisplacement', from {default}
  - **configurable as**: subfielddisplacement, displacement
* `pressure`: Pressure subfield.
  - **current value**: 'subfieldpressure', from {default}
  - **configurable as**: subfieldpressure, pressure

## Example

Example of setting `SolnDispPres` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem]
solution = pylith.problems.SolnDispPres
:::

