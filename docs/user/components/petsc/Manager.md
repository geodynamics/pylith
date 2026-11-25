% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# Manager

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.petsc.Manager`
:Journal name: `petsc`

Manage PETSc options.

## Example

Example of setting `Manager` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.petsc]
ts_monitor = true
ksp_monitor = true
ksp_converged_reason = true
snes_monitor = true
snes_converged_reason = true
snes_linesearch_monitor = true
:::

% End of file
