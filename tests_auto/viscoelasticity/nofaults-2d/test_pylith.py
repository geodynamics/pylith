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

from pylith.tests.FullTestApp import TestDriver, TestCase

import unittest
# import pdb


class TestApp(TestDriver):
    """
    Driver application for full-scale tests.
    """

    def __init__(self):
        """
        Constructor.
        """
        # pdb.set_trace()
        TestDriver.__init__(self)
        return

    def _suite(self):
        """
        Create test suite.
        """
        # pdb.set_trace()
        suite = unittest.TestSuite()

        import TestAxialTractionMaxwell
        for test in TestAxialTractionMaxwell.test_cases():
            suite.addTest(unittest.makeSuite(test))

        return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
    TestCase.parse_args()
    TestApp().main()


# End of file
