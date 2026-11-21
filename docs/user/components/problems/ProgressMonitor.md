% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# ProgressMonitor

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.ProgressMonitor`
:Journal name: `progressmonitor`

Abstract base class for simulation progress monitor.

## Pyre Properties

* `filename`=\<str\>: Name of output file.
  - **default value**: ''
  - **current value**: '', from {default}
* `update_percent`=\<float\>: Frequency of progress updates (percent).
  - **default value**: 5.0
  - **current value**: 5.0, from {default}
  - **validator**: (greater than 0)

% End of file
