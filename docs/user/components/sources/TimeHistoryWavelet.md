# TimeHistoryWavelet

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.sources.TimeHistoryWavelet`
:Journal name: `timehistorywavelet`

User-specified time history source time function.

This source time function allows specification of an arbitrary wavelet shape using a time history database.
The time history should contain normalized amplitudes (typically between 0 and 1) as a function of time.

Implements `SourceTimeFunctionMomentTensorForce`.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for source time function parameters.
  - **current value**: 'auxiliary_subfields', from {file='...', line=..., function='__set__'}
  - **configurable as**: auxiliary_subfields
* `time_history`: Time history with normalized amplitude as a function of time.
  - **current value**: 'nullcomponent', from {default}
  - **configurable as**: nullcomponent, time_history

## Pyre Properties

* `use_time_history`=\<bool\>: Use time history term in time-dependent expression.
  - **default value**: True
  - **current value**: True, from {default}

## Example

Example of setting `TimeHistoryWavelet` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.sources.source]
source_time_function = pylith.sources.TimeHistoryWavelet

# Specify time history database
source_time_function.time_history = spatialdata.spatialdb.TimeHistory
source_time_function.time_history.description = Source time history
source_time_function.time_history.filename = wavelet.timedb

# Source parameters specified in the auxiliary field database
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Source properties
db_auxiliary_field.values = [moment_tensor_xx, moment_tensor_yy, moment_tensor_xy, moment_tensor_zz, time_history_start_time]
db_auxiliary_field.data = [1.0e12*Pa*s, 1.0e12*Pa*s, 0.0*Pa*s, 1.0e12*Pa*s, 0.0*s]
:::

