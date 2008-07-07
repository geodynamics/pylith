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
class Integrator:

  def __init__(self, dt):
    self.dt = dt


  def stableTimeStep(self):
    return self.dt


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
    integrators = [Integrator(2.0*second),
                   Integrator(0.5*second)]

    # Set time step
    dt = 0.5*second / 1.2
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Increase stable time step, but time step should not change (skipped)
    integrators[1].dt = 0.8*second
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Reduce time step even if though should have skipped
    integrators[1].dt = 0.2*second
    dt = 0.2*second / 1.2
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Skip adjusting time step
    integrators[1].dt = 0.8*second
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Skip adjusting time step
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Adjust time step and stability factor
    tstep.stabilityFactor = 2.0
    dt = 0.8*second / 2.0
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Skip adjusting time step
    integrators[1].dt = 2.0*second
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Skip adjusting time step
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Adjust time step with value bigger than max
    dt = 1.0*second
    self.assertEqual(dt, tstep.timeStep(integrators))

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

    
    integrators = [Integrator(3.0*second),
                   Integrator(2.4*second)]
    dt = 2.4*second / 1.2
    tstep.timeStep(integrators)
    self.assertEqual(dt, tstep.currentStep())

    return


# End of file 
