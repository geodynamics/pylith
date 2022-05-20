# ProgressMonitorStep

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.ProgressMonitorStep`
:Journal name: `progressmonitorstep`

Progress monitor for problems with a given number of steps, such as Green's functions problem.

Implementes `ProgressMonitor`.

## Pyre Properties

* `filename`=\<str\>: Name of output file.
  - **default value**: 'progress.txt'
  - **current value**: 'progress.txt', from {default}
* `update_percent`=\<float\>: Frequency of progress updates (percent).
  - **default value**: 5.0
  - **current value**: 5.0, from {default}
  - **validator**: (greater than 0)

## Example

Example of setting `ProgressMonitorStep` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.timedependent.progress_monitor]
filename = output/greensfns01-progress.txt
:::

