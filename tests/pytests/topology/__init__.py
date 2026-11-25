# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


from . import (
    TestDistributor,
    TestField,
    TestMesh,
    TestRefineUniform,
    TestSubfield,
)


def test_modules():
    modules = [
        TestDistributor,
        TestField,
        TestMesh,
        TestRefineUniform,
        TestSubfield,
    ]
    return modules


# End of file
