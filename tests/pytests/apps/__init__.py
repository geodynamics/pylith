# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


from . import TestPyLithApp, TestEqInfoApp


def test_modules():
    return [
        TestPyLithApp,
        TestEqInfoApp,
    ]


# End of file
