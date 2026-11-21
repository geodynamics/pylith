% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# Physics

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.Physics`
:Journal name: `physics`

Abstract base class for objects defining physics.

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

% End of file
