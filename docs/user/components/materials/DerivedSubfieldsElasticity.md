# DerivedSubfieldsElasticity

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.DerivedSubfieldsElasticity`
:Journal name: `derivedsubfieldselasticity`

Derived subfields associated with the elasticity equation.

For elastic materials these derived subfields are available for output in addition to the solution and auxiliary subfields.

## Pyre Facilities

* `cauchy_strain`: Cauchy strain subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, cauchy_strain
* `cauchy_stress`: Cauchy stress subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, cauchy_stress

## Example

Example of setting `DerivedSubfieldsElasticity` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# The basis order for stress and strain should be at least 1 less than the basis order for displacement.
[pylithapp.problem.materials.mat_elastic.derived_subfields]
cauchy_stress.basis_order = 0
cauchy_strain.basis_order = 0
:::

