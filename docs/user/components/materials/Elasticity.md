# Elasticity

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.Elasticity`
:Journal name: `elasticity`

Material behavior governed by the elasticity equation.

Implements `Material`.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for auxiliary subfields.
  - **current value**: 'auxiliary_subfields', from {file='/Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pythia/pyre/inventory/ConfigurableClass.py', line=26, function='__set__'}
  - **configurable as**: auxiliary_subfields
* `bulk_rheology`: Bulk rheology for elastic material.
  - **current value**: 'isotropiclinearelasticity', from {default}
  - **configurable as**: isotropiclinearelasticity, bulk_rheology
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
* `use_body_force`=\<bool\>: Include body force term in elasticity equation.
  - **default value**: False
  - **current value**: False, from {default}

## Example

Example of setting `Elasticity` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.materials.mat_elastic]
label_value = 4
description = Upper crust elastic material
use_body_force = False
bulk_rheology = pylith.materials.IsotropicLinearElasticity

auxiliary_subfields.density.basis_order = 0
derived_subfields.cauchy_stress.basis_order = 1
derived_subfields.cauchy_strain.basis_order = 1
:::

