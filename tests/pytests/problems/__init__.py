# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
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
