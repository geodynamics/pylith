# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import unittest

from pylith.testing.TestCases import TestComponent, TestAbstractComponent, make_suite
import pylith.problems.Solution
import pylith.problems.SolnDisp
import pylith.problems.SolnDispLagrange
import pylith.problems.SolnDispPres
import pylith.problems.SolnDispPresLagrange
import pylith.problems.SolnDispPresTracStrain
import pylith.problems.SolnDispVel
import pylith.problems.SolnDispVelLagrange

class TestSolution(TestComponent):
    """Unit testing of Solution object.
    """
    _class = pylith.problems.Solution.Solution
    _factory = pylith.problems.Solution.solution


class TestSolnDisp(TestAbstractComponent):
    """Unit testing of SolnDisp object.
    """
    _class = pylith.problems.SolnDisp.SolnDisp


class TestSolnDispLagrange(TestAbstractComponent):
    """Unit testing of SolnDispLagrange object.
    """
    _class = pylith.problems.SolnDispLagrange.SolnDispLagrange


class TestSolutionDispLagrange(TestAbstractComponent):
    """Unit testing of Solution object.
    """
    _class = pylith.problems.SolnDispLagrange.Solution
    _factory = pylith.problems.SolnDispLagrange.solution


class TestSolnDispPres(TestAbstractComponent):
    """Unit testing of SolnDispPres object.
    """
    _class = pylith.problems.SolnDispPres.SolnDispPres


class TestSolutionDispPres(TestAbstractComponent):
    """Unit testing of Solution object.
    """
    _class = pylith.problems.SolnDispPres.Solution
    _factory = pylith.problems.SolnDispPres.solution


class TestSolnDispPresLagrange(TestAbstractComponent):
    """Unit testing of SolnDispPresLagrange object.
    """
    _class = pylith.problems.SolnDispPresLagrange.SolnDispPresLagrange


class TestSolutionDispPresLagrange(TestAbstractComponent):
    """Unit testing of Solution object.
    """
    _class = pylith.problems.SolnDispPresLagrange.Solution
    _factory = pylith.problems.SolnDispPresLagrange.solution


class TestSolnDispPresTracStrain(TestAbstractComponent):
    """Unit testing of SolnDispPresTracStrain object.
    """
    _class = pylith.problems.SolnDispPresTracStrain.SolnDispPresTracStrain


class TestSolutionDispPresTracStrain(TestAbstractComponent):
    """Unit testing of Solution object.
    """
    _class = pylith.problems.SolnDispPresTracStrain.Solution
    _factory = pylith.problems.SolnDispPresTracStrain.solution


class TestSolnDispVel(TestAbstractComponent):
    """Unit testing of SolnDispVel object.
    """
    _class = pylith.problems.SolnDispVel.SolnDispVel


class TestSolutionDispVel(TestAbstractComponent):
    """Unit testing of Solution object.
    """
    _class = pylith.problems.SolnDispVel.Solution
    _factory = pylith.problems.SolnDispVel.solution


class TestSolnDispVelLagrange(TestAbstractComponent):
    """Unit testing of SolnDispVelLagrange object.
    """
    _class = pylith.problems.SolnDispVelLagrange.SolnDispVelLagrange


class TestSolutionDispVelLagrange(TestAbstractComponent):
    """Unit testing of Solution object.
    """
    _class = pylith.problems.SolnDispVelLagrange.Solution
    _factory = pylith.problems.SolnDispVelLagrange.solution


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [
        TestSolution,
        TestSolnDisp,
        TestSolnDispLagrange,
        TestSolutionDispLagrange,
        TestSolnDispPres,
        TestSolutionDispPres,
        TestSolnDispPresLagrange,
        TestSolutionDispPresLagrange,
        TestSolnDispPresTracStrain,
        TestSolutionDispPresTracStrain,
        TestSolnDispVel,
        TestSolutionDispVel,
        TestSolnDispVelLagrange,
        TestSolutionDispVelLagrange,
    ]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
