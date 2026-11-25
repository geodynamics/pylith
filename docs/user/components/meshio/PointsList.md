% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# PointsList

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.PointsList`
:Journal name: `pointslist`

Reader for a list of points from an ASCII file.

:::{seealso}
See [`OutputSolnPoints` Component](OutputSolnPoints.md).
:::

## Pyre Facilities

* `coordsys`: Coordinate system associated with points.
  - **current value**: 'cscart', from {default}
  - **configurable as**: cscart, coordsys

## Pyre Properties

* `comment_delimiter`=\<str\>: Delimiter for comments.
  - **default value**: '#'
  - **current value**: '#', from {default}
* `filename`=\<str\>: Filename for list of points.
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateFilename at 0x11f3209d0>
* `value_delimiter`=\<str\>: Delimiter used to separate values.
  - **default value**: None
  - **current value**: None, from {default}

## Example

Example of setting `PointsList` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[points]
filename = stations.txt
comment_delimiter = #
value_delimiter = ,

coordsys = spatialdata.geocoords.CSCart
coordsys.space_dim = 2
:::

% End of file
