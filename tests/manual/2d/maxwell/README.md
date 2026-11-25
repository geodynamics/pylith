% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

This is a simple test of the MaxwellPlaneStrain viscoelastic model. The
results may be compared against the results from a 3D model with the
z-displacements constrained to zero.  There are four models:

Plane strain axial displacement problem:
pylith planestrain2d.cfg axialdisp2d.cfg

Plane strain shear displacement problem:
pylith planestrain2d.cfg sheardisp2d.cfg

3D simulated plane strain axial displacement problem:
pylith planestrain3d.cfg axialdisp3d.cfg

3D simulated plane strain shear displacement problem:
pylith planestrain3d.cfg sheardisp3d.cfg

The 2D and 3D problems should yield the same displacements, strains,
stresses, and viscous strains when the appropriate components are compared.

% End of file
