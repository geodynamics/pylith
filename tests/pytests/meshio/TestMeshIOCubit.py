#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/pytests/meshio/TestMeshIOCubit.py
#
# @brief Unit testing of Python MeshIOCubit object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.meshio.MeshIOCubit import (MeshIOCubit, mesh_io)


class TestMeshIOCubit(TestComponent):
    """Unit testing of MeshIOCubit object.
    """
    _class = MeshIOCubit
    _factory = mesh_io


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMeshIOCubit))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
