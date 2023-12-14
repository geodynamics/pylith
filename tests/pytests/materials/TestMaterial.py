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
# @file tests/pytests/materials/TestMaterial.py
#
# @brief Unit testing of Python TestMaterial object.

import unittest

from pylith.testing.UnitTestApp import TestAbstractComponent
from pylith.materials.Material import Material


class TestMaterial(TestAbstractComponent):
    """Unit testing of Material object.
    """
    _class = Material


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMaterial))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
