#!/usr/bin/env nemesis
# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


import numpy

from pylith.testing.FullTestApp import TestDriver
from TestEqInfo import TestEqInfo, run_eqinfo


class TestEqInfoTri(TestEqInfo):
    """Test suite for testing pylith_eqinfo with tri3 meshes."""

    def setUp(self):
        """Setup for test."""
        run_eqinfo("tri", ["tri.cfg"])
        return

    def test_stats(self):
        """Check fault stats."""
        import stats_tri

        timestamp = numpy.array([0.0, 1.0], dtype=numpy.float64)

        oneE = stats_tri.RuptureStats()
        oneE.timestamp = timestamp
        oneE.ruparea = numpy.array([1.5 + 2.0, 1.5 + 2.0], dtype=numpy.float64)
        slip0 = (0.2**2 + 0.5**2) ** 0.5
        slip1 = (0.5**2 + 0.4**2) ** 0.5
        oneE.potency = numpy.array([slip0 * 1.5 + slip1 * 2.0, 0.1 * 1.5 + 0.2 * 2.0], dtype=numpy.float64)
        oneE.moment = numpy.array(
            [slip0 * 1.5 * 1.0e10 + slip1 * 2.0 * 2.0e10, 0.1 * 1.5 * 1.0e10 + 0.2 * 2.0 * 2.0e10], dtype=numpy.float64
        )
        self._check(oneE, stats_tri.one)

        twoE = stats_tri.RuptureStats()
        twoE.timestamp = timestamp
        twoE.ruparea = numpy.array([1.5, 0.0], dtype=numpy.float64)
        twoE.potency = numpy.array([0.1 * 1.5, 0.0], dtype=numpy.float64)
        twoE.moment = numpy.array([0.1 * 1.5 * 1.0e10, 0.0], dtype=numpy.float64)
        self._check(twoE, stats_tri.two)

        all_E = stats_tri.RuptureStats()
        all_E.timestamp = timestamp
        all_E.ruparea = oneE.ruparea + twoE.ruparea
        all_E.potency = oneE.potency + twoE.potency
        all_E.moment = oneE.moment + twoE.moment
        self._check(all_E, stats_tri.all)
        return


# # -------------------------------------------------------------------------------------------------
def test_classes():
    return [
        TestEqInfoTri,
    ]


# ----------------------------------------------------------------------
if __name__ == "__main__":
    TestDriver(test_classes=test_classes()).run()


# End of file
