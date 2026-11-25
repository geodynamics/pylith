% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# Conventions

:::{warning}
This is a warning.
:::

:::{important}
This is something important.
:::

:::{tip}
This is a tip, helpful hint, or suggestion.
:::

For features recently added to PyLith, we show the version number when they
were added.
*New in v3.0.0.*

## Command Line Arguments

Example of a command line argument: `--help`.

## Filenames and Directories

Example of filenames and directories: `pylith`, `/usr/local`.

## Unix Shell Commands

Commands entered into a Unix shell (i.e., terminal) are shown in a box.
Comments are delimited by the # character. We use `$` to indicate the bash shell prompt.

```{code-block} bash
#This is a comment.
$ ls -l
```

## Excerpts of cfg Files

Example of an excerpt from a `.cfg` file:

```{code-block} cfg
# This is a comment.
[pylithapp.problem]
timestep = 2.0*s
bc = [x_pos, x_neg]
```

% End of file
