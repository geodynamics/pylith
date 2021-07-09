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
# @file tests/pytests/materials/TestPoroelasticity.py
#
# @brief Unit testing of Python TestPoroelasticity object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.materials.Poroelasticity import (Poroelasticity, material)


class TestPoroelasticity(TestComponent):
    """Unit testing of Poroelasticity object.
    """
    _class = Poroelasticity
    _factory = material


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestPoroelasticity))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
