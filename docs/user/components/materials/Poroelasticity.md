# Poroelasticity

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.Poroelasticity`
:Journal name: `poroelasticity`

Material behavior governed by the poroelasticity equation.

Implements `Material`.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for auxiliary subfields.
  - **current value**: 'auxiliary_subfields', from {file='/Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pythia/pyre/inventory/ConfigurableClass.py', line=26, function='__set__'}
  - **configurable as**: auxiliary_subfields
* `bulk_rheology`: Bulk rheology for poroelastic material.
  - **current value**: 'isotropiclinearporoelasticity', from {default}
  - **configurable as**: isotropiclinearporoelasticity, bulk_rheology
* `db_auxiliary_field`: Database for physical property parameters.
  - **current value**: 'simpledb', from {default}
  - **configurable as**: simpledb, db_auxiliary_field
* `derived_subfields`: Discretization of derived subfields.
  - **current value**: 'derived_subfields', from {file='/Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pythia/pyre/inventory/ConfigurableClass.py', line=26, function='__set__'}
  - **configurable as**: derived_subfields
* `observers`: Observers (e.g., output).
  - **current value**: 'singlephysicsobserver', from {default}
  - **configurable as**: singlephysicsobserver, observers

## Pyre Properties

* `description`=\<str\>: Descriptive label for material.
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateDescription at 0x1248ef790>
* `label`=\<str\>: Name of label for material. Currently only 'material-id' is allowed.
  - **default value**: 'material-id'
  - **current value**: 'material-id', from {default}
  - **validator**: (in ['material-id'])
* `label_value`=\<int\>: Value of label for material.
  - **default value**: 1
  - **current value**: 1, from {default}
* `use_body_force`=\<bool\>: Include body force term in Poroelasticity equation.
  - **default value**: False
  - **current value**: False, from {default}
* `use_source_density`=\<bool\>: Include source_density term in Poroelasticity equation.
  - **default value**: False
  - **current value**: False, from {default}
* `use_state_variables`=\<bool\>: Update porosity state variable using compaction formulation.
  - **default value**: False
  - **current value**: False, from {default}

## Example

Example of setting `Poroelasticity` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.materials.mat_poroelastic]
description = Upper crust poroelastic material
label_value = 3
use_body_force = True
use_source_density = False
use_state_variables = True
bulk_rheology = pylith.materials.IsotropicLinearPoroelasticity

auxiliary_subfields.density.basis_order = 0
auxiliary_subfields.body_force.basis_order = 0
:::

