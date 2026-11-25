% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# Python unit tests

We use the standard Python unittest module to implement the Python unit tests.

In reimplementing the Python unit tests in `tests/pytests`, we have setup bare bones tests that simply check that the objects are instantiated successfully, configure executes without errors, and factory functions return correct objects.
These tests are implemented in the `pylith.testing.UnitTestApp` Python object.
Tests of individual classes need only specify the class and the factory function.
We plan to expand these tests to include verification that information is properly transferred from Python to the underlying C++ object (if it exists).

% End of file
