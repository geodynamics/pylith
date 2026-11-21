% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# ConvertMeshApp

Full name: `pylith.apps.ConvertMeshApp`

Python ConvertMeshApp application.

## Pyre Facilities

* `reader`=\<component name\>: Reader for input mesh file.
  - **current value**: 'meshiopetsc', from {default}
  - **configurable as:** meshiopetsc, reader
* `writer`=\<component name\>: Writer for output mesh file.
  - **current value**: 'meshiopetsc', from {default}
  - **configurable as**: meshiopetsc, writer

## Pyre Properties

* `check_topology`=\<bool\>: Check topology of imported mesh.
  - **default value**: True
  - **current value**: True, from {default}

% End of file
