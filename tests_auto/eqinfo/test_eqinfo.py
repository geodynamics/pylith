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

import unittest

from pylith.tests.FullTestApp import TestDriver


class TestApp(TestDriver):
    """
    Driver application for full-scale tests.
    """

    def __init__(self):
        """
        Constructor.
        """
        TestDriver.__init__(self)
        return

    def _suite(self):
        """
        Create test suite.
        """
        suite = unittest.TestSuite()

        from TestEqInfoLine import TestEqInfoLine
        suite.addTest(unittest.makeSuite(TestEqInfoLine))

        from TestEqInfoTri import TestEqInfoTri
        suite.addTest(unittest.makeSuite(TestEqInfoTri))

        from TestEqInfoQuad import TestEqInfoQuad
        suite.addTest(unittest.makeSuite(TestEqInfoQuad))

        return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
    TestApp().main()


# End of file
