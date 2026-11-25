% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

This is a very simple test to make sure the stress/strain output for
HDF5 is consistent with output for VTK files. We use the same problem
setup for both cases (simple combination of Dirichlet and Neumann BC),
and then run the problem once using VTK output (output_vtk.cfg) and once
using HDF5 output (output_hdf5.cfg). The tensor components should match,
and should be consistent with what is in the PyLith manual.

% End of file
