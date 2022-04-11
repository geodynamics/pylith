# Material

% WARNING: Do not edit; this is a generated file!
Full name: `pylith.materials.Material`

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

* `id`=\<int\>: Material identifier (from mesh generator).
  - **default value**: 0
  - **current value**: 0, from {default}
* `label`=\<str\>: Descriptive label for material.
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateLabel at 0x117b40310>

