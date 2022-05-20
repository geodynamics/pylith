# KinSrcConstRate

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.faults.KinSrcConstRate`
:Journal name: `kinsrcconstrate`

Constant slip rate slip time function.

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

Example of setting `KinSrcConstRate` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
origin_time = 100*year

db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Constant slip rate slip time function auxiliary field spatial database
db_auxiliary_field.values = [initiation_time, slip_rate_left_lateral, slip_rate_opening]
db_auxiliary_field.data = [0.0*s, -2.0*mm/year, 0.0*mm/year]
:::

