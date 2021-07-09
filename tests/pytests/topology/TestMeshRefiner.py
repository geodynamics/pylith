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
# @file tests/pytests/topology/TestMeshRefiner.py
#
# @brief Unit testing of Python MeshRefiner object.

import unittest

from pylith.testing.UnitTestApp import TestAbstractComponent
from pylith.topology.MeshRefiner import MeshRefiner


class TestMeshRefiner(TestAbstractComponent):
    """Unit testing of MeshRefiner object.
    """
    _class = MeshRefiner


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMeshRefiner))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
