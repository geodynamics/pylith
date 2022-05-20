# AuxSubfieldsFault

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.faults.AuxSubfieldsFault`
:Journal name: `auxsubfieldfault`

Auxiliary subfields associated with a fault.

## Pyre Facilities

* `slip`: Slip subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, slip

## Example

Example of setting `AuxSubfieldsFault` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# We set the basis order to represent linear variations in the slip subfield.
[pylithapp.problem.interfaces.fault.auxiliary_fields]
slip.basis_order = 1
:::

