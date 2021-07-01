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
# @file tests/fullscale/eqinfo/TestEqInfoTri.py
#
# @brief Test suite for testing pylith_eqinfo with tri fault meshes.


import numpy

from TestEqInfo import TestEqInfo, run_eqinfo


class TestEqInfoTri(TestEqInfo):
    """Test suite for testing pylith_eqinfo with tri3 meshes.
    """

    def setUp(self):
        """Setup for test.
        """
        run_eqinfo("tri", ["tri.cfg"])
        return

    def test_stats(self):
        """Check fault stats.
        """
        import stats_tri

        timestamp = numpy.array([0.0, 1.0], dtype=numpy.float64)

        oneE = stats_tri.RuptureStats()
        oneE.timestamp = timestamp
        oneE.ruparea = numpy.array([1.5 + 2.0, 1.5 + 2.0], dtype=numpy.float64)
        slip0 = (0.2**2 + 0.5**2)**0.5
        slip1 = (0.5**2 + 0.4**2)**0.5
        oneE.potency = numpy.array(
            [slip0 * 1.5 + slip1 * 2.0, 0.1 * 1.5 + 0.2 * 2.0], dtype=numpy.float64)
        oneE.moment = numpy.array([slip0 * 1.5 * 1.0e+10 + slip1 * 2.0 * 2.0e+10,
                                   0.1 * 1.5 * 1.0e+10 + 0.2 * 2.0 * 2.0e+10], dtype=numpy.float64)
        self._check(oneE, stats_tri.one)

        twoE = stats_tri.RuptureStats()
        twoE.timestamp = timestamp
        twoE.ruparea = numpy.array([1.5, 0.0], dtype=numpy.float64)
        twoE.potency = numpy.array([0.1 * 1.5, 0.0], dtype=numpy.float64)
        twoE.moment = numpy.array(
            [0.1 * 1.5 * 1.0e+10, 0.0], dtype=numpy.float64)
        self._check(twoE, stats_tri.two)

        allE = stats_tri.RuptureStats()
        allE.timestamp = timestamp
        allE.ruparea = oneE.ruparea + twoE.ruparea
        allE.potency = oneE.potency + twoE.potency
        allE.moment = oneE.moment + twoE.moment
        self._check(allE, stats_tri.all)
        return


# ----------------------------------------------------------------------
if __name__ == '__main__':
    import unittest

    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestEqInfoTri))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
