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
