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
# @file tests/pytests/materials/TestIncompressibleElasticity.py
#
# @brief Unit testing of Python TestIncompressibleElasticity object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.materials.IncompressibleElasticity import (IncompressibleElasticity, material)


class TestIncompressibleElasticity(TestComponent):
    """Unit testing of IncompressibleElasticity object.
    """
    _class = IncompressibleElasticity
    _factory = material


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestIncompressibleElasticity))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
