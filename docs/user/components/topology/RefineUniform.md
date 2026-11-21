% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# RefineUniform

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.topology.RefineUniform`
:Journal name: `refineuniform`

Uniform global mesh refinement in parallel.

Implements `MeshRefiner`.

## Pyre Properties

* `levels`=\<int\>: Number of refinement levels.
  - **default value**: 0
  - **current value**: 0, from {default}
  - **validator**: (greater than or equal to 0)

## Example

Example of setting `RefineUniform` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# Refine mesh twice to reduce size of cell edges by a factor of 4.
[pylithapp.mesh_generator.refiner]
levels = 0
:::

% End of file
