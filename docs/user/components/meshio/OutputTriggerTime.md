% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# OutputTriggerTime

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.OutputTriggerTime`
:Journal name: `outputtriggertime`

Define how often output is written in terms of elapsed simulation time.

:::{tip}
Due to floating point roundoff, it is usually a good idea to use a value that is a fraction of a time step smaller than the desired value.
:::

Implements `OutputTrigger`.

## Pyre Properties

* `elapsed_time`=\<dimensional\>: Elapsed time between writes.
  - **default value**: 0*s
  - **current value**: 0*s, from {default}

## Example

Example of setting `OutputTriggerTime` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[output_trigger]
elapsed_time = 0.9999*year
:::

% End of file
