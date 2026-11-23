#!/usr/bin/env nemesis
# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


# @file tests/fullscale/eqinfo/TestEqInfoLine.py
#
# @brief Test suite for testing pylith_eqinfo with 1-D fault meshes.

import numpy

from pylith.testing.FullTestApp import TestDriver
from TestEqInfo import TestEqInfo, run_eqinfo


class TestEqInfoLine(TestEqInfo):
    """Test suite for testing pylith_eqinfo with 1-D fault meshes."""

    def setUp(self):
        """Setup for test."""
        run_eqinfo("line", ["line.cfg"])
        return

    def test_stats(self):
        """Check fault stats."""
        import stats_line

        timestamp = numpy.array([0.0, 1.0], dtype=numpy.float64)

        oneE = stats_line.RuptureStats()
        oneE.timestamp = timestamp
        oneE.ruparea = numpy.array([2.5, 1.5], dtype=numpy.float64)
        oneE.potency = numpy.array([0.7 * 1.0 + 0.9 * 1.5, 0.4 * 1.5], dtype=numpy.float64)
        oneE.moment = oneE.potency * 1.0e10
        self._check(oneE, stats_line.one)

        twoE = stats_line.RuptureStats()
        twoE.timestamp = timestamp
        area0 = (1.5**2 + 1.0**2) ** 0.5
        area1 = (1.0**2 + 1.0**2) ** 0.5
        twoE.ruparea = numpy.array([area0 + area1, area0], dtype=numpy.float64)
        twoE.potency = numpy.array([0.9 * area0 + 0.7 * area1, 0.3 * area0], dtype=numpy.float64)
        twoE.moment = twoE.potency * 1.0e10
        self._check(twoE, stats_line.two)

        all_E = stats_line.RuptureStats()
        all_E.timestamp = timestamp
        all_E.ruparea = oneE.ruparea + twoE.ruparea
        all_E.potency = oneE.potency + twoE.potency
        all_E.moment = oneE.moment + twoE.moment
        self._check(all_E, stats_line.all)
        return


# -------------------------------------------------------------------------------------------------
def test_classes():
    return [
        TestEqInfoLine,
    ]


# ----------------------------------------------------------------------
if __name__ == "__main__":
    TestDriver(test_classes=test_classes()).run()


# End of file
