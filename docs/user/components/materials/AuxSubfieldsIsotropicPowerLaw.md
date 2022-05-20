# AuxSubfieldsIsotropicPowerLaw

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.AuxSubfieldsIsotropicPowerLaw`
:Journal name: `auxfieldsisotropicpowerlaw`

Auxiliary subfields associated with the isotropic linear power law viscoelastic bulk rheology.

## Pyre Facilities

* `bulk_modulus`: Bulk modulus subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, bulk_modulus
* `deviatoric_stress`: Deviatoric stress subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, deviatoric_stress
* `power_law_exponent`: Power-law exponent subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, power_law_exponent
* `power_law_reference_strain_rate`: Power-law reference strain rate subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, power_law_reference_strain_rate
* `power_law_reference_stress`: Power-law reference stress subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, power_law_reference_stress
* `reference_strain`: Reference strain subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, reference_strain
* `reference_stress`: Reference stress subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, reference_stress
* `shear_modulus`: Shear modulus subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, shear_modulus
* `viscous_strain`: Viscous strain subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, viscous_strain

## Example

Example of setting `AuxSubfieldsIsotropicPowerLaw` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.materials.mat_powerlaw.rheology.auxiliary_fields]
shear_modulus.basis_order = 1
bulk_modulus.basis_order = 1
power_law_reference_strain_rate = 1
power_law_reference_stress = 1
power_law_exponent.basis_order = 1
viscous_strain.basis_order = 1
stress.basis_order = 1
reference_stress.basis_order = 0
reference_strain.basis_order = 0
:::

