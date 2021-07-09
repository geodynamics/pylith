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
# @file tests/pytests/bc/TestDirichletTimeDependent.py
#
# @brief Unit testing of Python DirichletTimeDependent object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.bc.DirichletTimeDependent import (DirichletTimeDependent, boundary_condition)


class TestDirichletTimeDependent(TestComponent):
    """Unit testing of DirichletTimeDependent object.
    """
    _class = DirichletTimeDependent
    _factory = boundary_condition

    @staticmethod
    def customizeInventory(obj):
        obj.inventory.constrainedDOF = ["0", "2"]


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDirichletTimeDependent))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
