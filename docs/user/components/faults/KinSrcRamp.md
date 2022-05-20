# KinSrcRamp

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.faults.KinSrcRamp`
:Journal name: `kinsrcramp`

Linear ramp slip time function.

Implements `KinSrc`.

## Pyre Facilities

* `db_auxiliary_field`: Database for slip time function parameters.
  - **current value**: 'simpledb', from {default}
  - **configurable as**: simpledb, db_auxiliary_field

## Pyre Properties

* `origin_time`=\<dimensional\>: Origin time for slip source.
  - **default value**: 0*s
  - **current value**: 0*s, from {default}

## Example

Example of setting `KinSrcRamp` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
origin_time = 10*year

db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Ramp slip time function auxiliary field spatial database
db_auxiliary_field.values = [initiation_time, rise_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [0.0*s, 3.0*s, -2.0*m, 0.0*m]
:::

