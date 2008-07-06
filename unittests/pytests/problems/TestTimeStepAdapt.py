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

## @file unittests/pytests/problems/TestTimeStepAdapt.py

## @brief Unit testing of TimeStepAdapt object.

import unittest
from pylith.problems.TimeStepAdapt import TimeStepAdapt

from pyre.units.time import second

# ----------------------------------------------------------------------
class TestTimeStepAdapt(unittest.TestCase):
  """
  Unit testing of TimeStepAdapt object.
  """

  def test_numTimeSteps(self):
    """
    Test numTimeSteps().
    """
    tstep = TimeStepAdapt()
    tstep._configure()

    self.assertEqual(1, tstep.numTimeSteps())

    tstep.totalTime = 4.0*second
    tstep.maxDt = 2.0*second
    self.assertEqual(3, tstep.numTimeSteps())

    return


  def test_timeStep(self):
    """
    Test timeStep().
    """
    tstep = TimeStepAdapt()
    tstep._configure()
    
    tstep.adaptSkip = 2

    # Set time step
    dt = 0.5*second
    self.assertEqual(dt, tstep.timeStep(dt))

    # Time step should not change
    self.assertEqual(dt, tstep.timeStep(2.0*second))

    # Reduce time step even if should have skipped
    dt = 0.2*second
    self.assertEqual(dt, tstep.timeStep(dt))

    # Skip adjusting time step
    self.assertEqual(dt, tstep.timeStep(0.5*second))

    # Skip adjusting time step
    self.assertEqual(dt, tstep.timeStep(0.5*second))

    # Adjust time step
    dt = 0.8*second
    self.assertEqual(dt, tstep.timeStep(dt))

    # Skip adjusting time step
    self.assertEqual(dt, tstep.timeStep(2.0*second))

    # Skip adjusting time step
    self.assertEqual(dt, tstep.timeStep(2.0*second))

    # Adjust time step with value bigger than max
    dt = 1.0*second
    self.assertEqual(dt, tstep.timeStep(2.0*second))

    return


  def test_currentStep(self):
    """
    Test currentStep().
    """
    tstep = TimeStepAdapt()
    tstep._configure()
    tstep.maxDt = 10.0*second
    tstep.dt = tstep.maxDt
    
    self.assertEqual(10.0*second, tstep.currentStep())

    dt = 4.0*second
    tstep.timeStep(dt)
    self.assertEqual(dt, tstep.currentStep())

    return


# End of file 
