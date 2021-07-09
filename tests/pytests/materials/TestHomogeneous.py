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
# @file tests/pytests/materials/TestHomogeneous.py
#
# @brief Unit testing of Python Homogeneous object.

import unittest

from pylith.testing.UnitTestApp import TestAbstractComponent
from pylith.materials.Homogeneous import Homogeneous


class TestHomogeneous(TestAbstractComponent):
    """Unit testing of Homogeneous object.
    """
    _class = Homogeneous


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestHomogeneous))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
