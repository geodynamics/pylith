# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


from . import (
    TestMaterial,
    TestElasticity,
    TestIncompressibleElasticity,
    TestPoroelasticity,
    TestRheologies,
    TestAuxiliarySubfields,
    TestDerivedSubfields,
    TestHomogeneous,
)


def test_modules():
    return [
        TestMaterial,
        TestElasticity,
        TestIncompressibleElasticity,
        TestPoroelasticity,
        TestRheologies,
        TestAuxiliarySubfields,
        TestDerivedSubfields,
        TestHomogeneous,
    ]


# End of file
