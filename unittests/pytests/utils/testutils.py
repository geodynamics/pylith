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

# @file unittests/utils/testutils.py

# @brief Python application for testing utils code.

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

        from TestPylith import TestPylith
        suite.addTest(unittest.makeSuite(TestPylith))

        from TestEventLogger import TestEventLogger
        suite.addTest(unittest.makeSuite(TestEventLogger))

        from TestPylithVersion import TestPylithVersion
        suite.addTest(unittest.makeSuite(TestPylithVersion))

        from TestPetscVersion import TestPetscVersion
        suite.addTest(unittest.makeSuite(TestPetscVersion))

        from TestDependenciesVersion import TestDependenciesVersion
        suite.addTest(unittest.makeSuite(TestDependenciesVersion))

        from TestCollectVersionInfo import TestCollectVersionInfo
        suite.addTest(unittest.makeSuite(TestCollectVersionInfo))

        from TestConstants import TestConstants
        suite.addTest(unittest.makeSuite(TestConstants))

        return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
    app = TestApp()
    app.run()


# End of file
