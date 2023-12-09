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
