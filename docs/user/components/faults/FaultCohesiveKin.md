# FaultCohesiveKin

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.faults.FaultCohesiveKin`
:Journal name: `faultcohesivekin`

Fault surface with kinematic (prescribed) slip implemented with cohesive cells.

The fault may have an arbitrary number of kinematic sources for coseismic slip and creep.
They are superimposed at each time step to create the prescribed slip on the fault.

Implements `FaultCohesive`.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for auxiliary subfields.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, auxiliary_subfields
* `db_auxiliary_field`: (no documentation available)
  - **current value**: 'nullcomponent', from {default}
  - **configurable as**: nullcomponent, db_auxiliary_field
* `derived_subfields`: Discretization of derived subfields.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, derived_subfields
* `eq_ruptures`: Kinematic earthquake sources information.
  - **current value**: 'singlerupture', from {default}
  - **configurable as**: singlerupture, eq_ruptures
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

## Example

Example of setting `FaultCohesiveKin` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# Specify prescribed slip on a fault via two earthquakes in a 2D domain.
[pylithapp.problem.interfaces.fault]
label = fault
edge = fault_edge

observers.observer.data_fields = [slip]

# Two earthquakes with different slip time functions.
eq_ruptures = [quake10, quake50]
quake10 = pylith.faults.KinSrcBrune
quake50 = pylith.faults.KinSrcLiuCosine

# Rupture parameters for the first earthquake.
[pylithapp.problem.interfaces.fault.eq_ruptures.quake10]
origin_time = 10*year

db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Fault rupture auxiliary field spatial database
db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [0.0*s, -2.0*m, 0.0*m]

# Rupture parameters for the second earthquake.
[pylithapp.problem.interfaces.fault.eq_ruptures.quake50]
origin_time = 50*year

db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Fault rupture auxiliary field spatial database
db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [0.0*s, -1.0*m, 0.0*m]
:::

