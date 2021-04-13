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
# @file tests/pytests/problems/TestSolnDisp.py
#
# @brief Unit testing of Python SolnDisp object.

import unittest

from pylith.testing.UnitTestApp import TestAbstractComponent
from pylith.problems.SolnDisp import SolnDisp


class TestSolnDisp(TestAbstractComponent):
    """Unit testing of SolnDisp object.
    """
    _class = SolnDisp


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSolnDisp))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
