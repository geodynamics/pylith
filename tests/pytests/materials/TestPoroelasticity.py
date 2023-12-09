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
