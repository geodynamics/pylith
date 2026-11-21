# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


from . import TestCommunicator, TestReduce


def test_modules():
    modules = [
        TestCommunicator,
        TestReduce,
    ]
    return modules


# End of file
