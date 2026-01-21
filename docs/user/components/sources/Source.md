# Source

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.sources.Source`
:Journal name: `source`

Abstract base class for point sources.

This is the base class for all point source implementations in PyLith.
Point sources are used to specify internal forcing terms at discrete points within the domain.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for auxiliary subfields.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, auxiliary_subfields
* `db_auxiliary_field`: Database for source parameters.
  - **current value**: 'simpledb', from {default}
  - **configurable as**: simpledb, db_auxiliary_field
* `derived_subfields`: Discretization of derived subfields.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, derived_subfields
* `observers`: Observers (e.g., output).
  - **current value**: 'singlephysicsobserver', from {default}
  - **configurable as**: singlephysicsobserver, observers
* `reader`: Reader for source points list.
  - **current value**: 'pointslist', from {default}
  - **configurable as**: pointslist, reader

## Pyre Properties

* `description`=\<str\>: Descriptive label for source.
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateDescription at ...>
* `field`=\<str\>: Solution subfield associated with source.
  - **default value**: 'displacement'
  - **current value**: 'displacement', from {default}
* `label`=\<str\>: Name of label for source.
  - **default value**: 'source-id'
  - **current value**: 'source-id', from {default}
* `label_value`=\<int\>: Value of label identifying source.
  - **default value**: 1
  - **current value**: 1, from {default}

