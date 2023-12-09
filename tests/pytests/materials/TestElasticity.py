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
