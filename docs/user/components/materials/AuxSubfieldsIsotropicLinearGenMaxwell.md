# AuxSubfieldsIsotropicLinearGenMaxwell

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.AuxSubfieldsIsotropicLinearGenMaxwell`
:Journal name: `auxfieldsisotropiclineargenmaxwell`

Auxiliary subfields associated with the isotropic generalized Maxwell viscoelastic bulk rheology.

## Pyre Facilities

* `bulk_modulus`: Bulk modulus subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, bulk_modulus
* `maxwell_time`: Maxwell time subfield for 3 Maxwell elements.
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
* `shear_modulus_ratio`: Shear modulus ratio subfield for 3 Maxwell elements.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, shear_modulus_ratio
* `total_strain`: Total strain subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, total_strain
* `viscous_strain`: Viscous strain subfield for 3 Maxwell elements.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, viscous_strain

## Example

Example of setting `AuxSubfieldsIsotropicLinearGenMaxwell` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.materials.mat_genmaxwell.rheology.auxiliary_fields]
shear_modulus.basis_order = 1
bulk_modulus.basis_order = 1
maxwell_time.basis_order = 1
shear_modulus_ratio.basis_order = 1
total_strain.basis_order = 1
viscous_strain.basis_order = 1
reference_stress.basis_order = 0
reference_strain.basis_order = 0
:::

