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
