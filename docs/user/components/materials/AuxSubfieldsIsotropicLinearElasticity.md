# AuxSubfieldsIsotropicLinearElasticity

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.AuxSubfieldsIsotropicLinearElasticity`
:Journal name: `auxsubfieldsisotropiclinearelasticity`

Auxiliary subfields associated with the isotropic linear elastic bulk rheology.

:::{important}
The auxiliary subfields (internal representation of material properties) do not necessarily match the values in the spatial database.
For example, the spatial database uses density, Vp, and Vs instead of density, shear modulus, and bulk modulus because that is how they are usually characterized in seismic velocity models.
PyLith converts the values provided by the user in a spatial database to the internal representation stored in the auxiliary field.
:::

## Pyre Facilities

* `bulk_modulus`: Bulk modulus subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, bulk_modulus
* `reference_strain`: Reference strain subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, reference_strain
* `reference_stress`: Reference stress subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, reference_stress
* `shear_modulus`: Shear modulus subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, shear_modulus

## Example

Example of setting `AuxSubfieldsIsotropicLinearElasticity` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.materials.mat_elastic.rheology.auxiliary_fields]
shear_modulus.basis_order = 1
bulk_modulus.basis_order = 1
reference_stress.basis_order = 0
reference_strain.basis_order = 0
:::

