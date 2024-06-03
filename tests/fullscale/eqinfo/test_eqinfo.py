#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
import unittest

from pylith.testing.FullTestApp import TestDriver


class TestApp(TestDriver):
    """Driver application for full-scale tests.
    """

    def __init__(self):
        """Constructor.
        """
        TestDriver.__init__(self)
        return

    def _suite(self):
        """Create test suite.
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
