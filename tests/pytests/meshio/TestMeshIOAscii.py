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
# @file tests/pytests/meshio/TestMeshIOAscii.py
#
# @brief Unit testing of Python MeshIOAscii object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.meshio.MeshIOAscii import (MeshIOAscii, mesh_io)


class TestMeshIOAscii(TestComponent):
    """Unit testing of MeshIOAscii object.
    """
    _class = MeshIOAscii
    _factory = mesh_io


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMeshIOAscii))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
