% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# MeshReordering

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.initializers.MeshReordering`
:Journal name: `reorder_mesh`

Reorder mesh using the reverse Cuthill-McKee algorithm.

Implements `InitializePhase`.

## Example

Example of setting `MeshReordering` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.mesh_initializer.phases.reorder_mesh]
# No parameters
:::

% End of file
