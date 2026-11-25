% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

This test focuses on verifying that PyLith correctly handles spatial
variations in physical properties and updating state variables.

Fault slip is imposed at time t=0 and then the domain undergoes
viscoelastic relaxation with a short relaxation on one side of the
fault (completely relaxes) and a longer relaxation time on the other
side of the fault.

% End of file
