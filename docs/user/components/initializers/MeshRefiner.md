% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# MeshRefiner

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.initializers.MeshRefiner`
:Journal name: `mesh_refiner`

Refine mesh.

Implements `InitializePhase`.

## Pyre Facilities

* `refiner`: Mesh refiner.
  - **current value**: 'refineuniform', from {default}
  - **configurable as**: refineuniform, refiner

## Example

Example of setting `MeshRefiner` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.mesh_initializer.phases.refine_mesh]
refiner = pylith.topology.RefineUniform
:::

% End of file
