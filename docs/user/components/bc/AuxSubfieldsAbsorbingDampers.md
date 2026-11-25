% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# AuxSubfieldsAbsorbingDampers

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.bc.AuxSubfieldsAbsorbingDampers`
:Journal name: `auxsubfieldsabsorbingdampers`

Auxiliary subfields for the absorbing dampers boundary condition.

## Pyre Facilities

* `density`: Mass density subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, density
* `vp`: Dilatational (P) wave speed subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, vp
* `vs`: Shear (S) wave speed subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, vs

## Example

Example of setting `AuxSubfieldsAbsorbingDampers` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[absorbing_dampers_auxiliary_subfields]
density.basis_order = 0
vp.basis_order = 0
vs.basis_order = 0            
:::

% End of file
