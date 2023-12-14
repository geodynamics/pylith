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
# @file tests/pytests/topology/TestDistributor.py
#
# @brief Unit testing of Python Distributor object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.topology.Distributor import (Distributor, mesh_distributor)


class TestDistributor(TestComponent):
    """Unit testing of Distributor object.
    """
    _class = Distributor
    _factory = mesh_distributor


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDistributor))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
