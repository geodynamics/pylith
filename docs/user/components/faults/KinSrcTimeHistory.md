# KinSrcTimeHistory

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.faults.KinSrcTimeHistory`
:Journal name: `kinsrctimehistory`

Time history slip time function.

The slip time function is given in a `TimeHistory` database.

Implements `KinSrc`.

## Pyre Facilities

* `db_auxiliary_field`: Database for slip time function parameters.
  - **current value**: 'simpledb', from {default}
  - **configurable as**: simpledb, db_auxiliary_field
* `time_history`: Time history with normalized amplitude as a function of time.
  - **current value**: 'timehistory', from {default}
  - **configurable as**: timehistory, time_history

## Pyre Properties

* `origin_time`=\<dimensional\>: Origin time for slip source.
  - **default value**: 0*s
  - **current value**: 0*s, from {default}

## Example

Example of setting `KinSrcTimeHistory` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
origin_time = 10*year

db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Time history slip time function auxiliary field spatial database
db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [0.0*s, -2.0*m, 0.0*m]

time_history.description = Slip time function time history
time_history.filename = slipfn.timedb
:::

