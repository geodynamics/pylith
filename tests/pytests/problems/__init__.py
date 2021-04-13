from .TestInitialCondition import TestInitialCondition
from .TestInitialConditionDomain import TestInitialConditionDomain
from .TestInitialConditionPatch import TestInitialConditionPatch
from .TestPhysics import TestPhysics
from .TestProblem import TestProblem
from .TestTimeDependent import TestTimeDependent
from .TestProblemDefaults import TestProblemDefaults
from .TestProgressMonitor import TestProgressMonitor
from .TestProgressMonitorTime import TestProgressMonitorTime
from .TestSingleObserver import (TestSinglePhysicsObserver, TestSingleSolnObserver)
from .TestSolution import TestSolution
from .TestSolnDisp import TestSolnDisp
from .TestSolnDispLagrange import (TestSolnDispLagrange, TestSolutionDispLagrange)
from .TestSolnDispPres import (TestSolnDispPres, TestSolutionDispPres)
from .TestSolnDispLagrange import (TestSolnDispLagrange, TestSolutionDispLagrange)
from .TestSolnDispTracStrain import (TestSolnDispPresTracStrain, TestSolutionDispTracStrain)
from .TestSolnDispVel import (TestSolnDispVel, TestSolutionDispVel)
from .TestSolnDispVelLagrange import (TestSolnDispVelLagrange, TestSolutionDispVelLagrange)
from .TestSolutionSubfield import TestSolutionSubfield
from .TestSubfieldDisplacement import TestSubfieldDisplacement
from .TestSubfieldLagrangeFault import TestSubfieldLagrangeFault
from .TestSubfieldPressure import TestSubfieldPressure
from .TestSubfieldTemperature import TestSubfieldTemperature
from .TestSubfieldTraceStrain import TestSubfieldTraceStrain
from .TestSubfieldVelocity import TestSubfieldVelocity


def test_classes():
    classes = [
        TestInitialCondition,
        TestInitialConditionDomain,
        TestInitialConditionPatch,
        TestPhysics,
        TestProblem,
        TestTimeDependent,
        TestProblemDefaults,
        TestProgressMonitorTime,
        TestProgressMonitorTime,
        TestSingleSolnObserver,
        TestSinglePhysicsObserver,
        TestSolution,
        TestSolnDisp,
        TestSolnDispLagrange,
        TestSolutionDispLagrange,
        TestSolnDispPresTracStrain, 
        TestSolutionDispTracStrain,
        TestSolnDispVel, 
        TestSolutionDispVel,
        TestSolnDispVelLagrange, 
        TestSolutionDispVelLagrange,
        TestSolutionSubfield,
        TestSubfieldDisplacement,
        TestSubfieldLagrangeFault,
        TestSubfieldPressure,
        TestSubfieldTemperature,
        TestSubfieldTraceStrain,
        TestSubfieldVelocity,
    ]
    return classes


# End of file
