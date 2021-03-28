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
# @file tests/pytests/materials/TestElasticity.py
#
# @brief Unit testing of Python TestElasticity object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.materials.Elasticity import (Elasticity, material)


class TestElasticity(TestComponent):
    """Unit testing of Elasticity object.
    """
    _class = Elasticity
    _factory = material


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestElasticity))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
