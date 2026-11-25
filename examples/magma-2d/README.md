% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# Examples: 2D Magma Reservoir Using Poroelasticity

This example demonstrates the use of poroelasticity. It includes a magma reservoir with
poroelastic properties that differ from the surrounding domain.

## Meshing

We provide mesh files generated using Gmsh and Cubit. The parameter files are setup to
use the mesh from Gmsh.

## Step 1: Inflation

We impose flow into the conduit at the external boundary, which leads to inflation of the
magma reservoir.

To run the example:

```bash
pylith step01_inflation.cfg
```

## Step 2: Inflation with porosity state variable

Same as Step 1 with evolution of porosity state variable.

To run the example:

```bash
pylith step02_inflation_statevars.cfg
```

% End of file
