# SolnDispLagrange

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.SolnDispLagrange`
:Journal name: `solndisplagrange`

Container for solution subfields with displacement and fault Lagrange multiplier subfields.

## Pyre Facilities

* `displacement`: Displacement subfield.
  - **current value**: 'subfielddisplacement', from {default}
  - **configurable as**: subfielddisplacement, displacement
* `lagrange_fault`: Fault Lagrange multiplier subfield.
  - **current value**: 'subfieldlagrangefault', from {default}
  - **configurable as**: subfieldlagrangefault, lagrange_fault

## Example

Example of setting `SolnDispLagrange` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem]
solution = pylith.problems.SolnDispLagrange
:::

