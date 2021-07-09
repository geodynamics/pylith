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
