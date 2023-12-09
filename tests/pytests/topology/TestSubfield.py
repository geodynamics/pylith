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
# @file tests/pytests/topology/TestSubfield.py
#
# @brief Unit testing of Python Subfield object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.topology.Subfield import (Subfield, subfield)


class TestSubfield(TestComponent):
    """Unit testing of Subfield object.
    """
    _class = Subfield
    _factory = subfield


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSubfield))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
