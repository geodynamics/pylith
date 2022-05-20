# FaultCohesiveImpulses

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.faults.FaultCohesiveImpulses`
:Journal name: `faultcohesiveimpulses`

Fault surface with slip impulses for Green's functions implemented with cohesive cells.

The comopnents 

Implements `FaultCohesiveKin`.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for auxiliary subfields.
  - **current value**: 'auxiliary_subfields', from {file='/Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pythia/pyre/inventory/ConfigurableClass.py', line=26, function='__set__'}
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
* `impulse_dof`=\<list\>: Indices of impulse components; 0=fault opening, 1=left lateral, 2=reverse (3D only).
  - **default value**: []
  - **current value**: [], from {default}
  - **validator**: <function validateDOF at 0x124b029d0>
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
threshold = 0.5

# Create impulses at all points on the fault by specifying a uniform amplitude of 1.0.
# Impulses will be applied at any location with a slip component greater than the threshold.
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Slip impulse amplitude
db_auxiliary_field.values = [slip_left_lateral, slip_opening]
db_auxiliary_field.data = [1.0*m, 0.0*m]

# Represent the impulse as a linear variation in slip centered on each point.
auxiliary_subfields.slip.basis_order = 1
:::

