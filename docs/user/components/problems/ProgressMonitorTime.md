% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# ProgressMonitorTime

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.ProgressMonitorTime`
:Journal name: `progressmonitortime`

Progress monitor for time-dependent problem.

If you do not set the filename for the progress monitor, then PyLith will create one using the
simulation name from the application defaults settings.

## Pyre Properties

* `filename`=\<str\>: Name of output file.
  - **default value**: ''
  - **current value**: '', from {default}
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

% End of file
