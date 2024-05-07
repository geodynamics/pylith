# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import unittest

from pylith.testing.TestCases import (TestAbstractComponent, TestComponent, make_suite)
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


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [
        TestKinSrc,
        TestKinSrcConstRate,
        TestKinSrcStep,
        TestKinSrcRamp,
        TestKinSrcBrune,
        TestKinSrcLiuCos,
        TestKinSrcTimeHistory,
    ]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
