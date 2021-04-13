#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#
# @file tests/pytests/bc/TestNeumannTimeDependent.py
#
# @brief Unit testing of Python NeumannTimeDependent object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.bc.NeumannTimeDependent import (NeumannTimeDependent, boundary_condition)


class TestNeumannTimeDependent(TestComponent):
    """Unit testing of NeumannTimeDependent object.
    """
    _class = NeumannTimeDependent
    _factory = boundary_condition


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestNeumannTimeDependent))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
