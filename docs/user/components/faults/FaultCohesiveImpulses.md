# FaultCohesiveImpulses

% WARNING: Do not edit; this is a generated file!
Full name: `pylith.faults.FaultCohesiveImpulses`

Fault surface with slip impulses for Green's functions implemented with cohesive cells.

Implements `FaultCohesiveKin`.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for auxiliary subfields.
  - **current value**: 'auxiliary_subfields', from {file='/software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pythia/pyre/inventory/ConfigurableClass.py', line=26, function='__set__'}
  - **configurable as**: auxiliary_subfields
* `db_auxiliary_field`: (no documentation available)
  - **current value**: 'nullcomponent', from {default}
  - **configurable as**: nullcomponent, db_auxiliary_field
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
* `impulse_dof`=\<list\>: Indices of impulse components (0=1st DOF, 1=2nd DOF, etc).
  - **default value**: []
  - **current value**: [], from {default}
  - **validator**: <function validateDOF at 0x7fe78c8fe820>
* `label`=\<str\>: Name of label identifier for fault.
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateLabel at 0x7fe78c8fea60>
* `label_value`=\<int\>: Value of label identifier for fault.
  - **default value**: 1
  - **current value**: 1, from {default}
* `ref_dir_1`=\<list\>: First choice for reference direction to discriminate among tangential directions in 3-D.
  - **default value**: [0.0, 0.0, 1.0]
  - **current value**: [0.0, 0.0, 1.0], from {default}
  - **validator**: <function validateDir at 0x7fe78c8feaf0>
* `ref_dir_2`=\<list\>: Second choice for reference direction to discriminate among tangential directions in 3-D.
  - **default value**: [0.0, 1.0, 0.0]
  - **current value**: [0.0, 1.0, 0.0], from {default}
  - **validator**: <function validateDir at 0x7fe78c8feaf0>
* `threshold`=\<dimensional\>: Threshold for non-zero amplitude.
  - **default value**: 1e-06*m
  - **current value**: 1e-06*m, from {default}
  - **validator**: (greater than or equal to 0*m)

## Example

Example of setting `FaultCohesiveImpulses` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.greensfns]
interfaces = [fault]
interfaces.fault = pylith.faults.FaultCohesiveImpulses

[pylithapp.greensfns.interfaces.fault]
label = fault
label_value = 20

# Impulses for left-lateral slip (dof=1)
impulse_dof = [1]
:::

