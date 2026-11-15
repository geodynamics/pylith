# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
from . import (
    TestInitializer,
    TestSerial,
    TestParallel,
    TestConvert,
    TestInitializePhase,
    TestMeshReader,
    TestMeshWriter,
    TestMeshReordering,
    TestMeshRefiner,
    TestMeshDistributor,
    TestMeshInsertInterfaces,
)


def test_modules():
    return [
        TestInitializer,
        TestSerial,
        TestParallel,
        TestConvert,
        TestInitializePhase,
        TestMeshReader,
        TestMeshWriter,
        TestMeshReordering,
        TestMeshRefiner,
        TestMeshDistributor,
        TestMeshInsertInterfaces,
    ]


# End of file
