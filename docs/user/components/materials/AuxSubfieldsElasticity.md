# AuxSubfieldsElasticity

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.AuxSubfieldsElasticity`
:Journal name: `auxsubfieldselasticity`

Auxiliary subfields associated with the elasticity equation.

Setting the parameters for a subfield does not turn on its use.
The [`Elasticity` Component](Elasticity.md) has flags for including or excluding terms in the elasticity equation.

## Pyre Facilities

* `body_force`: Body force subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, body_force
* `density`: Density subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, density
* `gravitational_acceleration`: Gravitational acceleration subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, gravitational_acceleration

## Example

Example of setting `AuxSubfieldsElasticity` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# We set the basis order to represent linear variations in the density and body 
# force subfields and a uniform gravitational acceleration subfield.
[pylithapp.problem.materials.mat_elastic.auxiliary_fields]
density.basis_order = 1
body_force.basis_order = 1
gravitational_acceleration.basis_order = 0
:::

