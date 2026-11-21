% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# SolnDispPresTracStrain

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.SolnDispPresTracStrain`
:Journal name: `solndispprestracstrain`

Container for solution subfields with displacement, pore pressure, and trace strain subfields.

## Pyre Facilities

* `displacement`: Displacement subfield.
  - **current value**: 'subfielddisplacement', from {default}
  - **configurable as**: subfielddisplacement, displacement
* `pressure`: Pressure subfield.
  - **current value**: 'subfieldpressure', from {default}
  - **configurable as**: subfieldpressure, pressure
* `trace_strain`: Trace strain subfield.
  - **current value**: 'subfieldtracestrain', from {default}
  - **configurable as**: subfieldtracestrain, trace_strain

## Example

Example of setting `SolnDispPresTracStrain` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem]
solution = pylith.problems.SolnDispPresTracStrain
:::

% End of file
