% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# MeshReader

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.initializers.MeshReader`
:Journal name: `mesh_reader`

Read mesh.

Implements `InitializePhase`.

## Pyre Facilities

* `reader`: Mesh reader.
  - **current value**: 'meshiopetsc', from {default}
  - **configurable as**: meshiopetsc, reader

## Example

Example of setting `MeshReader` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.mesh_initializer.phases.read_mesh]
reader = pylith.meshio.MeshIOPetsc
:::

% End of file
