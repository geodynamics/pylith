% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# BoundaryCondition

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.bc.BoundaryCondition`
:Journal name: `boundarycondition`

Abstract base class for boundary conditions.

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

* `field`=\<str\>: Solution subfield associated with boundary condition.
  - **default value**: 'displacement'
  - **current value**: 'displacement', from {default}
* `label`=\<str\>: Name of label identifying boundary.
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateLabel at 0x124bbc4c0>
* `label_value`=\<int\>: Value of label identifying boundary (tag of physical group in Gmsh files).
  - **default value**: 1
  - **current value**: 1, from {default}

% End of file
