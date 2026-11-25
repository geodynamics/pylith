% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# DataWriterHDF5Ext

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.DataWriterHDF5Ext`
:Journal name: `datawriterhdf5`

Writer of solution, auxiliary, and derived subfields to an HDF5 file with datasets stored in external binary files.

Implements `DataWriter`.

## Pyre Properties

* `filename`=\<str\>: Name of HDF5 file.
  - **default value**: ''
  - **current value**: '', from {default}

## Example

Example of setting `DataWriterHDF5Ext` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[data_writer]
filename = domain_solution.h5
:::

% End of file
