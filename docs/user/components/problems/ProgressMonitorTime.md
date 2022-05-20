# ProgressMonitorTime

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.ProgressMonitorTime`
:Journal name: `progressmonitortime`

Progress monitor for time-dependent problem.

## Pyre Properties

* `filename`=\<str\>: Name of output file.
  - **default value**: 'progress.txt'
  - **current value**: 'progress.txt', from {default}
* `t_units`=\<str\>: Units used for simulation time in output.
  - **default value**: 'year'
  - **current value**: 'year', from {default}
* `update_percent`=\<float\>: Frequency of progress updates (percent).
  - **default value**: 5.0
  - **current value**: 5.0, from {default}
  - **validator**: (greater than 0)

## Example

Example of setting `ProgressMonitorTime` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.timedependent.progress_monitor]
filename = output/step01-progress.txt
t_units = year
:::

