# AuxSubfieldsMomentTensorForce

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.sources.AuxSubfieldsMomentTensorForce`
:Journal name: `auxsubfieldsmomenttensorforce`

Container for moment tensor force auxiliary subfields.

This component defines the discretization for the moment tensor and time delay auxiliary subfields.

## Pyre Facilities

* `moment_tensor`: Moment tensor subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, moment_tensor
* `time_delay`: Time delay subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, time_delay

## Example

Example of setting `AuxSubfieldsMomentTensorForce` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.sources.source]
auxiliary_subfields.moment_tensor.basis_order = 0
auxiliary_subfields.time_delay.basis_order = 0
:::

