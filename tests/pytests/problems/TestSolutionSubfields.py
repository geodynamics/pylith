# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import unittest

from pylith.testing.TestCases import TestComponent, make_suite

import pylith.problems.SolutionSubfield
import pylith.problems.SubfieldDisplacement
import pylith.problems.SubfieldLagrangeFault
import pylith.problems.SubfieldPressure
import pylith.problems.SubfieldTemperature
import pylith.problems.SubfieldTraceStrain
import pylith.problems.SubfieldVelocity

class TestSolutionSubfield(TestComponent):
    """Unit testing of SolutionSubfield object.
    """
    _class = pylith.problems.SolutionSubfield.SolutionSubfield
    _factory = pylith.problems.SolutionSubfield.soln_subfield


class TestSubfieldDisplacement(TestComponent):
    """Unit testing of SubfieldDisplacement object.
    """
    _class = pylith.problems.SubfieldDisplacement.SubfieldDisplacement
    _factory = pylith.problems.SubfieldDisplacement.soln_subfield


class TestSubfieldLagrangeFault(TestComponent):
    """Unit testing of SubfieldLagrangeFault object.
    """
    _class = pylith.problems.SubfieldLagrangeFault.SubfieldLagrangeFault
    _factory = pylith.problems.SubfieldLagrangeFault.soln_subfield


class TestSubfieldPressure(TestComponent):
    """Unit testing of SubfieldPressure object.
    """
    _class = pylith.problems.SubfieldPressure.SubfieldPressure
    _factory = pylith.problems.SubfieldPressure.soln_subfield


class TestSubfieldTemperature(TestComponent):
    """Unit testing of SubfieldTemperature object.
    """
    _class = pylith.problems.SubfieldTemperature.SubfieldTemperature
    _factory = pylith.problems.SubfieldTemperature.soln_subfield


class TestSubfieldTraceStrain(TestComponent):
    """Unit testing of SubfieldTraceStrain object.
    """
    _class = pylith.problems.SubfieldTraceStrain.SubfieldTraceStrain
    _factory = pylith.problems.SubfieldTraceStrain.soln_subfield


class TestSubfieldVelocity(TestComponent):
    """Unit testing of SubfieldVelocity object.
    """
    _class = pylith.problems.SubfieldVelocity.SubfieldVelocity
    _factory = pylith.problems.SubfieldVelocity.soln_subfield


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [
        TestSolutionSubfield,
        TestSubfieldDisplacement,
        TestSubfieldLagrangeFault,
        TestSubfieldPressure,
        TestSubfieldTemperature,
        TestSubfieldTraceStrain,
        TestSubfieldVelocity,
    ]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
