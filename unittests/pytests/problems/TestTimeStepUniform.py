#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

## @file unittests/pytests/problems/TestTimeStepUniform.py

## @brief Unit testing of TimeStepUniform object.

import unittest
from pylith.problems.TimeStepUniform import TimeStepUniform

from pyre.units.time import second

# ----------------------------------------------------------------------
class TestTimeStepUniform(unittest.TestCase):
  """
  Unit testing of TimeStepUniform object.
  """

  def test_numTimeSteps(self):
    """
    Test numTimeSteps().
    """
    tstep = TimeStepUniform()
    tstep._configure()

    self.assertEqual(1, tstep.numTimeSteps())

    tstep.totalTime = 4.0*second
    tstep.dt = 2.0*second
    self.assertEqual(3, tstep.numTimeSteps())

    return


  def test_timeStep(self):
    """
    Test timeStep().
    """
    tstep = TimeStepUniform()
    tstep._configure()

    integrators = None
    mesh = None

    self.assertEqual(1.0*second, tstep.timeStep(mesh, integrators))

    tstep.dt = 1.0e-4*second
    self.assertEqual(1.0e-4*second, tstep.timeStep(mesh, integrators))

    return


  def test_currentStep(self):
    """
    Test currentStep().
    """
    tstep = TimeStepUniform()
    tstep._configure()

    self.assertEqual(1.0*second, tstep.currentStep())

    tstep.dt = 1.0e-4*second
    self.assertEqual(1.0e-4*second, tstep.currentStep())

    return


# End of file 
