#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
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

import unittest

from pylith.testing.UnitTestApp import TestAbstractComponent
from pylith.utils.CollectVersionInfo import CollectVersionInfo


class TestCollectVersionInfo(TestAbstractComponent):
    _class = CollectVersionInfo

    def test_asDict(self):
        info = CollectVersionInfo.asDict()
        self.assertTrue(info)
        return

    def test_asString(self):
        s = CollectVersionInfo.asString()
        self.assertTrue(s)
        return

if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestCollectVersionInfo))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
