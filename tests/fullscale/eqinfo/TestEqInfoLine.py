#!/usr/bin/env nemesis
#
# ----------------------------------------------------------------------
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
# ----------------------------------------------------------------------
#
# @file tests/fullscale/eqinfo/TestEqInfoLine.py
#
# @brief Test suite for testing pylith_eqinfo with 1-D fault meshes.

import numpy

from TestEqInfo import TestEqInfo, run_eqinfo


class TestEqInfoLine(TestEqInfo):
    """Test suite for testing pylith_eqinfo with 1-D fault meshes.
    """

    def setUp(self):
        """Setup for test.
        """
        run_eqinfo("line", ["line.cfg"])
        return

    def test_stats(self):
        """Check fault stats.
        """
        import stats_line

        timestamp = numpy.array([0.0, 1.0], dtype=numpy.float64)

        oneE = stats_line.RuptureStats()
        oneE.timestamp = timestamp
        oneE.ruparea = numpy.array([2.5, 1.5], dtype=numpy.float64)
        oneE.potency = numpy.array(
            [0.7 * 1.0 + 0.9 * 1.5, 0.4 * 1.5], dtype=numpy.float64)
        oneE.moment = oneE.potency * 1.0e+10
        self._check(oneE, stats_line.one)

        twoE = stats_line.RuptureStats()
        twoE.timestamp = timestamp
        area0 = (1.5**2 + 1.0**2)**0.5
        area1 = (1.0**2 + 1.0**2)**0.5
        twoE.ruparea = numpy.array([area0 + area1, area0], dtype=numpy.float64)
        twoE.potency = numpy.array(
            [0.9 * area0 + 0.7 * area1, 0.3 * area0], dtype=numpy.float64)
        twoE.moment = twoE.potency * 1.0e+10
        self._check(twoE, stats_line.two)

        allE = stats_line.RuptureStats()
        allE.timestamp = timestamp
        allE.ruparea = oneE.ruparea + twoE.ruparea
        allE.potency = oneE.potency + twoE.potency
        allE.moment = oneE.moment + twoE.moment
        self._check(allE, stats_line.all)
        return


# ----------------------------------------------------------------------
if __name__ == '__main__':
    import unittest

    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestEqInfoLine))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
