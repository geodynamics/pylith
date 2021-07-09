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
