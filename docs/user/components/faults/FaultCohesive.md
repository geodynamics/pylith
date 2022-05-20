# FaultCohesive

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.faults.FaultCohesive`
:Journal name: `fault`

Abstract base class for a fault surface implemeted with cohesive cells.

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

* `edge`=\<str\>: Name of label identifier for buried fault edges.
  - **default value**: ''
  - **current value**: '', from {default}
* `edge_value`=\<int\>: Value of label identifier for buried fault edges.
  - **default value**: 1
  - **current value**: 1, from {default}
* `label`=\<str\>: Name of label identifier for fault.
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateLabel at 0x124b02ca0>
* `label_value`=\<int\>: Value of label identifier for fault.
  - **default value**: 1
  - **current value**: 1, from {default}
* `ref_dir_1`=\<list\>: First choice for reference direction to discriminate among tangential directions in 3-D.
  - **default value**: [0.0, 0.0, 1.0]
  - **current value**: [0.0, 0.0, 1.0], from {default}
  - **validator**: <function validateDir at 0x124b02dc0>
* `ref_dir_2`=\<list\>: Second choice for reference direction to discriminate among tangential directions in 3-D.
  - **default value**: [0.0, 1.0, 0.0]
  - **current value**: [0.0, 1.0, 0.0], from {default}
  - **validator**: <function validateDir at 0x124b02dc0>

