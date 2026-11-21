# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


from . import (
    TestDirichletTimeDependent,
    TestNeumannTimeDependent,
    TestAbsorbingDampers,
    TestAuxiliarySubfields,
    TestZeroDB,
)


def test_modules():
    return [
        TestDirichletTimeDependent,
        TestNeumannTimeDependent,
        TestAbsorbingDampers,
        TestAuxiliarySubfields,
        TestZeroDB,
    ]


# End of file
