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
# @file tests/pytests/problems/TestProgressMonitorTime.py
#
# @brief Unit testing of Python ProgressMonitorTime object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.problems.ProgressMonitorTime import (ProgressMonitorTime, progress_monitor)


class TestProgressMonitorTime(TestComponent):
    """Unit testing of ProgressMonitorTime object.
    """
    _class = ProgressMonitorTime
    _factory = progress_monitor


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestProgressMonitorTime))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
