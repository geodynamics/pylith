# AuxSubfieldsIsotropicLinearMaxwell

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.AuxSubfieldsIsotropicLinearMaxwell`
:Journal name: `auxfieldsisotropiclinearmaxwell`

Auxiliary subfields associated with the isotropic linear Maxwell viscoelastic bulk rheology.

## Pyre Facilities

* `bulk_modulus`: Bulk modulus subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, bulk_modulus
* `maxwell_time`: Maxwell time subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, maxwell_time
* `reference_strain`: Reference strain subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, reference_strain
* `reference_stress`: Reference stress subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, reference_stress
* `shear_modulus`: Shear modulus subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, shear_modulus
* `total_strain`: Total strain subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, total_strain
* `viscous_strain`: Viscous strain subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, viscous_strain

## Example

Example of setting `AuxSubfieldsIsotropicLinearMaxwell` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.materials.mat_maxwell.rheology.auxiliary_fields]
shear_modulus.basis_order = 1
bulk_modulus.basis_order = 1
maxwell_time.basis_order = 1
total_strain.basis_order = 1
viscous_strain.basis_order = 1
reference_stress.basis_order = 0
reference_strain.basis_order = 0
:::

