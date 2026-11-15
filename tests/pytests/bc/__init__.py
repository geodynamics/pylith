# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
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
