# Material

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.Material`
:Journal name: `material`

Abstract base class for a bulk material.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for auxiliary subfields.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, auxiliary_subfields
* `db_auxiliary_field`: Database for physical property parameters.
  - **current value**: 'simpledb', from {default}
  - **configurable as**: simpledb, db_auxiliary_field
* `derived_subfields`: Discretization of derived subfields.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, derived_subfields
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

