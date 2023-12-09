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
