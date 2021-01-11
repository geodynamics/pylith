#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

# @file tests/topology/testdriver.py

# @brief Python application for testing topology code.

from pylith.tests.UnitTestApp import UnitTestApp

import unittest


class TestApp(UnitTestApp):
    """
    Test application.
    """

    def __init__(self):
        """
        Constructor.
        """
        UnitTestApp.__init__(self)
        return

    def _suite(self):
        """
        Setup the test suite.
        """

        suite = unittest.TestSuite()

        from TestMesh import TestMesh
        suite.addTest(unittest.makeSuite(TestMesh))

        from TestSubmesh import TestSubmesh
        suite.addTest(unittest.makeSuite(TestSubmesh))

        from TestFieldBase import TestFieldBase
        suite.addTest(unittest.makeSuite(TestFieldBase))

        from TestMeshField import TestMeshField
        suite.addTest(unittest.makeSuite(TestMeshField))

        from TestMeshFields import TestMeshFields
        suite.addTest(unittest.makeSuite(TestMeshFields))

        from TestSolutionFields import TestSolutionFields
        suite.addTest(unittest.makeSuite(TestSolutionFields))

        from TestJacobian import TestJacobian
        suite.addTest(unittest.makeSuite(TestJacobian))

        from TestMeshGenerator import TestMeshGenerator
        suite.addTest(unittest.makeSuite(TestMeshGenerator))

        from TestMeshImporter import TestMeshImporter
        suite.addTest(unittest.makeSuite(TestMeshImporter))

        from TestRefineUniform import TestRefineUniform
        suite.addTest(unittest.makeSuite(TestRefineUniform))

        return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
    app = TestApp()
    app.run()


# End of file
