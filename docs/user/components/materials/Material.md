% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# Material

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.Material`
:Journal name: `material`

Abstract base class for a bulk material.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for auxiliary subfields.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, auxiliary_subfields
* `db_auxiliary_field`: Database for physical property parameters.
  - **current value**: 'simpledb', from {default}
  - **configurable as**: simpledb, db_auxiliary_field
* `derived_subfields`: Discretization of derived subfields.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, derived_subfields
* `observers`: Observers (e.g., output).
  - **current value**: 'singlephysicsobserver', from {default}
  - **configurable as**: singlephysicsobserver, observers

## Pyre Properties

* `description`=\<str\>: Descriptive label for material (deprecated in v5.0; will be removed in v6.0).
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateDescription at 0x10608f2e0>
* `label`=\<str\>: Name of label for material. Currently only 'material-id' is allowed.
  - **default value**: 'material-id'
  - **current value**: 'material-id', from {default}
  - **validator**: (in ['material-id'])
* `label_value`=\<int\>: Value of label for material.
  - **default value**: 1
  - **current value**: 1, from {default}

% End of file
