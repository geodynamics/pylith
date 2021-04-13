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
