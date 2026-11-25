% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# AuxSubfieldsFault

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.faults.AuxSubfieldsFault`
:Journal name: `auxsubfieldfault`

Auxiliary subfields associated with a fault.

## Pyre Facilities

* `slip`: Slip subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, slip

## Example

Example of setting `AuxSubfieldsFault` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# We set the basis order to represent linear variations in the slip subfield.
[pylithapp.problem.interfaces.fault.auxiliary_fields]
slip.basis_order = 1
:::

% End of file
