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
# @file tests/pytests/faults/TestKinSrc.py
#
# @brief Unit testing of Python KinSrc object.

import unittest

from pylith.testing.UnitTestApp import (TestAbstractComponent, TestComponent)
import pylith.faults.KinSrc
import pylith.faults.KinSrcConstRate
import pylith.faults.KinSrcStep
import pylith.faults.KinSrcRamp
import pylith.faults.KinSrcBrune
import pylith.faults.KinSrcLiuCos
import pylith.faults.KinSrcTimeHistory


class TestKinSrc(TestAbstractComponent):
    """Unit testing of KinSrc object.
    """
    _class = pylith.faults.KinSrc.KinSrc


class TestKinSrcConstRate(TestComponent):
    """Unit testing of KinSrcConstRate object.
    """
    _class = pylith.faults.KinSrcConstRate.KinSrcConstRate
    _factory = pylith.faults.KinSrcConstRate.eq_kinematic_src


class TestKinSrcStep(TestComponent):
    """Unit testing of KinSrcStep object.
    """
    _class = pylith.faults.KinSrcStep.KinSrcStep
    _factory = pylith.faults.KinSrcStep.eq_kinematic_src


class TestKinSrcRamp(TestComponent):
    """Unit testing of KinSrcRamp object.
    """
    _class = pylith.faults.KinSrcRamp.KinSrcRamp
    _factory = pylith.faults.KinSrcRamp.eq_kinematic_src


class TestKinSrcBrune(TestComponent):
    """Unit testing of KinSrcBrune object.
    """
    _class = pylith.faults.KinSrcBrune.KinSrcBrune
    _factory = pylith.faults.KinSrcBrune.eq_kinematic_src


class TestKinSrcLiuCos(TestComponent):
    """Unit testing of KinSrcLiuCos object.
    """
    _class = pylith.faults.KinSrcLiuCos.KinSrcLiuCos
    _factory = pylith.faults.KinSrcLiuCos.eq_kinematic_src


class TestKinSrcTimeHistory(TestComponent):
    """Unit testing of KinSrcTimeHistory object.
    """
    _class = pylith.faults.KinSrcTimeHistory.KinSrcTimeHistory
    _factory = pylith.faults.KinSrcTimeHistory.eq_kinematic_src


if __name__ == "__main__":
    suite = unittest.TestSuite()
    for cls in [
        TestKinSrc,
        TestKinSrcConstRate,
        TestKinSrcStep,
        TestKinSrcRamp,
        TestKinSrcBrune,
        TestKinSrcLiuCos,
        TestKinSrcTimeHistory,
    ]:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
