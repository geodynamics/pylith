% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# SolnDispPresVel

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.SolnDispPresVel`
:Journal name: `solndisppresvel`

Container for solution subfields with displacement, pore pressure, and velocity subfields.

## Pyre Facilities

* `displacement`: Displacement subfield.
  - **current value**: 'subfielddisplacement', from {default}
  - **configurable as**: subfielddisplacement, displacement
* `pressure`: Pressure subfield.
  - **current value**: 'subfieldpressure', from {default}
  - **configurable as**: subfieldpressure, pressure
* `velocity`: Velocity subfield.
  - **current value**: 'subfieldvelocity', from {default}
  - **configurable as**: subfieldvelocity, velocity

## Example

Example of setting `SolnDispPresVel` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem]
solution = pylith.problems.SolnDispPresVel
:::

% End of file
