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
# @file tests/pytests/meshio/TestOutputPhysics.py
#
# @brief Unit testing of Python OutputPhysics object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.meshio.OutputPhysics import (OutputPhysics, observer)


class TestOutputPhysics(TestComponent):
    """Unit testing of OutputPhysics object.
    """
    _class = OutputPhysics
    _factory = observer


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestOutputPhysics))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
