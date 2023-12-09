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
# @file tests/pytests/topology/TestMeshImporter.py
#
# @brief Unit testing of Python MeshImporter object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.topology.MeshImporter import (MeshImporter, mesh_generator)


class TestMeshImporter(TestComponent):
    """Unit testing of MeshImporter object.
    """
    _class = MeshImporter
    _factory = mesh_generator


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMeshImporter))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
