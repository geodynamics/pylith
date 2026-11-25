% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# DerivedSubfieldsPoroelasticity

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.DerivedSubfieldsPoroelasticity`
:Journal name: `derivedsubfieldsporoelasticity`

Derived subfields associated with the poroelasticity equation.

For poroelastic materials these derived subfields are available for output in addition to the solution and auxiliary subfields.

## Pyre Facilities

* `bulk_density`: Bulk density subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, bulk_density
* `cauchy_strain`: Cauchy strain subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, cauchy_strain
* `cauchy_stress`: Cauchy stress subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, cauchy_stress
* `water_content`: Water content subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, water_content

## Example

Example of setting `DerivedSubfieldsPoroelasticity` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# The basis order for stress and strain should be at least 1 less than the basis order for displacement.
[pylithapp.problem.materials.mat_elastic.derived_subfields]
cauchy_stress.basis_order = 0
cauchy_strain.basis_order = 0
bulk_density.basis_order = 1
water_content.basis_order = 1
:::

% End of file
