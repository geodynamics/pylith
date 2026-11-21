% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# Step 1: Gravitational Body Forces

This example involves a static simulation that solves for the deformation from loading by gravitational body forces.
{numref}`fig:example:troubleshooting:2d:step01:diagram` shows the boundary conditions on the domain.
The parameters specific to this example are in `step01_gravity.cfg`.

:::{figure-md} fig:example:troubleshooting:2d:step01:diagram
<img src="figs/step01-diagram.*" alt="" scale="75%">

We apply roller boundary conditions on the lateral sides and bottom of the domain.
:::

See [examples/reverse-2d Step 1](../reverse-2d/step01-gravity.md) for more information.

% End of file
