# AuxSubfieldsPoroelasticity

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.AuxSubfieldsPoroelasticity`
:Journal name: `auxsubfieldsporoelasticity`

Auxiliary subfields associated with the poroelasticity equation.

Setting the parameters for a subfield does not turn on its use.
The [`Poroelasticity` Component](Poroelasticity.md) has flags for including or excluding terms in the poroelasticity equation.

## Pyre Facilities

* `body_force`: Body force subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, body_force
* `fluid_density`: Fluid density subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, fluid_density
* `fluid_viscosity`: Fluid viscosity subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, fluid_viscosity
* `gravitational_acceleration`: Gravitational acceleration subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, gravitational_acceleration
* `porosity`: Porosity subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, porosity
* `solid_density`: Solid density subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, solid_density
* `source_density`: Source density subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, source_density

## Example

Example of setting `AuxSubfieldsPoroelasticity` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.materials.mat_poroelastic.auxiliary_fields]
porosity.basis_order = 1
solid_density.basis_order = 1
fluid_density.basis_order = 0
fluid_viscosity.basis_order = 1
body_force.basis_order = 1
source_density.basis_order = 1
gravitational_acceleration.basis_order = 0
:::

