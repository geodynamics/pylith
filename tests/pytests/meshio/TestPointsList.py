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
# @file tests/pytests/meshio/TestPointsList.py
#
# @brief Unit testing of Python PointsList object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.meshio.PointsList import (PointsList, points_list)


class TestPointsList(TestComponent):
    """Unit testing of PointsList object.
    """
    _class = PointsList
    _factory = points_list


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestPointsList))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
