% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# AuxSubfieldsIsotropicLinearPoroelasticity

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.AuxSubfieldsIsotropicLinearPoroelasticity`
:Journal name: `auxfieldsisotropiclinearporoelasticity`

Auxiliary subfields associated with the isotropic linear poroelastic bulk rheology.

## Pyre Facilities

* `biot_coefficient`: Biot coefficient subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, biot_coefficient
* `biot_modulus`: Biot modulus subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, biot_modulus
* `drained_bulk_modulus`: Drained bulk modulus subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, drained_bulk_modulus
* `isotropic_permeability`: Isotropic permeability subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, isotropic_permeability
* `reference_strain`: Reference strain subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, reference_strain
* `reference_stress`: Reference stress subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, reference_stress
* `shear_modulus`: Shear modulus subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, shear_modulus
* `tensor_permeability`: Tensor permeability subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, tensor_permeability

## Example

Example of setting `AuxSubfieldsIsotropicLinearPoroelasticity` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.materials.mat_poroelastic.rheology.auxiliary_fields]
shear_modulus.basis_order = 1
biot_coefficient.basis_order = 0
isotropic_permeability.basis_order = 0
drained_bulk_modulus.basis_order = 1
biot_modulus.basis_order = 1
:::

% End of file
