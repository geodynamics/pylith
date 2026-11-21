% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# MeshWriter

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.initializers.MeshWriter`
:Journal name: `mesh_writer`

Write mesh.

Implements `InitializePhase`.

## Pyre Facilities

* `writer`: Mesh writer.
  - **current value**: 'meshiopetsc', from {default}
  - **configurable as**: meshiopetsc, writer

## Example

Example of setting `MeshWriter` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.mesh_initializer.phases.write_mesh]
writer = pylith.meshio.MeshIOPetsc
:::

% End of file
