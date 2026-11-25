% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# Convert

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.initializers.Convert`
:Journal name: `convert_phases`

Mesh initialization phases for reordering and converting to a different format.

:::{seealso}
See [`Initializer` Component](Initializer.md).
:::

## Pyre Facilities

* `read_mesh`: Read mesh.
  - **current value**: 'mesh_reader', from {default}
  - **configurable as**: mesh_reader, read_mesh
* `reorder_mesh`: Reorder mesh using reverse Cuthill-McKee algorithm.
  - **current value**: 'reorder_mesh', from {default}
  - **configurable as**: reorder_mesh
* `write_mesh`: Write mesh.
  - **current value**: 'mesh_writer', from {default}
  - **configurable as**: mesh_writer, write_mesh

% End of file
