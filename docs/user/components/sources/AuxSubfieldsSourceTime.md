# AuxSubfieldsSourceTime

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.sources.AuxSubfieldsSourceTime`
:Journal name: `auxsubfieldssourcetime`

Container for source time function auxiliary subfields.

This component defines the discretization for the center frequency auxiliary subfield used by wavelet source time functions.

## Pyre Facilities

* `center_frequency`: Center frequency subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, center_frequency

## Example

Example of setting `AuxSubfieldsSourceTime` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.sources.source]
source_time_function.auxiliary_subfields.center_frequency.basis_order = 0
:::

