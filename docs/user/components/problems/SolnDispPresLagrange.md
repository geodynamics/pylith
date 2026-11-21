% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# SolnDispPresLagrange

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.SolnDispPresLagrange`
:Journal name: `solndisppres`

Container for solution subfields with displacement, pressure, and fault Lagrange multiplier subfields.

## Pyre Facilities

* `displacement`: Displacement subfield.
  - **current value**: 'subfielddisplacement', from {default}
  - **configurable as**: subfielddisplacement, displacement
* `lagrange_multiplier_fault`: Fault Lagrange multiplier subfield.
  - **current value**: 'subfieldlagrangefault', from {default}
  - **configurable as**: subfieldlagrangefault, lagrange_multiplier_fault
* `pressure`: Pressure subfield.
  - **current value**: 'subfieldpressure', from {default}
  - **configurable as**: subfieldpressure, pressure

## Example

Example of setting `SolnDispPresLagrange` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem]
solution = pylith.problems.SolnDispPresLagrange
:::

% End of file
