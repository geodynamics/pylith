#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
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

from pylith.utils.CollectVersionInfo import CollectVersionInfo


class TestCollectVersionInfo(unittest.TestCase):

    def test_asDict(self):
        info = CollectVersionInfo.asDict()
        self.assertTrue(info)
        return

    def test_asString(self):
        s = CollectVersionInfo.asString()
        self.assertTrue(s)
        return

# End of file
