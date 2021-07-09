#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ======================================================================
#
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
