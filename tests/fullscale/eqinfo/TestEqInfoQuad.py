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
# @file tests/fullscale/eqinfo/TestEqInfoQuad.py
#
# @brief Test suite for testing pylith_eqinfo with quad fault meshes.

import numpy

from TestEqInfo import TestEqInfo, run_eqinfo


class TestEqInfoQuad(TestEqInfo):
    """Test suite for testing pylith_eqinfo with quad4 meshes.
    """

    def setUp(self):
        """Setup for test.
        """
        run_eqinfo("quad", ["quad.cfg"])
        return

    def test_stats(self):
        """Check fault stats.
        """
        import stats_quad

        timestamp = numpy.array([5.0], dtype=numpy.float64)

        area0 = 1.5 * 1.75
        area1 = 1.35 * 1.5
        slip0 = (1.0**2 + 1.2**2)**0.5
        slip1 = (1.4**2 + 1.6**2)**0.5
        oneE = stats_quad.RuptureStats()
        oneE.timestamp = timestamp
        oneE.ruparea = numpy.array([area0 + area1], dtype=numpy.float64)
        oneE.potency = numpy.array(
            [slip0 * area0 + slip1 * area1], dtype=numpy.float64)
        oneE.moment = numpy.array(
            [slip0 * area0 * 1.0e+10 + area1 * slip1 * 2.0e+10], dtype=numpy.float64)
        self._check(oneE, stats_quad.one)

        area0 = 1.5 * 1.75
        area1 = 1.35 * 1.5
        slip0 = (1.0**2 + 1.2**2)**0.5
        slip1 = (1.4**2 + 1.6**2)**0.5
        twoE = stats_quad.RuptureStats()
        twoE.timestamp = timestamp
        twoE.ruparea = numpy.array([area0 + area1], dtype=numpy.float64)
        twoE.potency = numpy.array(
            [slip0 * area0 + slip1 * area1], dtype=numpy.float64)
        twoE.moment = numpy.array(
            [slip0 * area0 * 1.0e+10 + area1 * slip1 * 2.0e+10], dtype=numpy.float64)
        self._check(twoE, stats_quad.two)

        allE = stats_quad.RuptureStats()
        allE.timestamp = timestamp
        allE.ruparea = oneE.ruparea + twoE.ruparea
        allE.potency = oneE.potency + twoE.potency
        allE.moment = oneE.moment + twoE.moment
        self._check(allE, stats_quad.all)
        return


# ----------------------------------------------------------------------
if __name__ == '__main__':
    import unittest

    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestEqInfoQuad))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
