# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


from . import (
    TestInitialCondition,
    TestInitialConditionDomain,
    TestInitialConditionPatch,
    TestPhysics,
    TestProblem,
    TestTimeDependent,
    TestProblemDefaults,
    TestProgressMonitor,
    TestProgressMonitorTime,
    TestSingleObserver,
    TestSolution,
    TestSolutionSubfields,
)


def test_modules():
    modules = [
        TestInitialCondition,
        TestInitialConditionDomain,
        TestInitialConditionPatch,
        TestPhysics,
        TestProblem,
        TestTimeDependent,
        TestProblemDefaults,
        TestProgressMonitor,
        TestProgressMonitorTime,
        TestSingleObserver,
        TestSolution,
        TestSolutionSubfields,
    ]
    return modules


# End of file
