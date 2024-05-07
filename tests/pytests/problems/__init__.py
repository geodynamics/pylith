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
