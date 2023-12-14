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
# @file tests/pytests/topology/TestRefineUniform.py
#
# @brief Unit testing of Python RefineUniform object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.topology.RefineUniform import (RefineUniform, mesh_refiner)


class TestRefineUniform(TestComponent):
    """Unit testing of RefineUniform object.
    """
    _class = RefineUniform
    _factory = mesh_refiner


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestRefineUniform))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
