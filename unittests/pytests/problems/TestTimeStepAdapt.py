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

  def setUp(self):
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()
    normalizer.setTimeScale(2.0*second)

    tstep = TimeStepAdapt()
    tstep._configure()
    tstep.preinitialize()
    tstep.initialize(normalizer)
    self.tstep = tstep
    return
  

  def test_initialize(self):
    """
    Test initialize().
    """
    tstep = self.tstep

    self.assertEqual(0.0, tstep.totalTime)
    self.assertEqual(0.5, tstep.maxDt)
    return
  

  def test_numTimeSteps(self):
    """
    Test numTimeSteps().
    """
    tstep = self.tstep

    self.assertEqual(1, tstep.numTimeSteps())

    tstep.totalTime = 4.0
    tstep.maxDt = 2.0
    self.assertEqual(3, tstep.numTimeSteps())

    return


  def test_timeStep(self):
    """
    Test timeStep().
    """
    tstep = self.tstep

    tstep.adaptSkip = 2
    integrators = [Integrator(2.0),
                   Integrator(0.5)]

    # Set time step
    dt = 0.5 / 1.2
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Increase stable time step, but time step should not change (skipped)
    integrators[1].dt = 0.8
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Reduce time step even if though should have skipped
    integrators[1].dt = 0.2
    dt = 0.2 / 1.2
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Skip adjusting time step
    integrators[1].dt = 0.8
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Skip adjusting time step
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Adjust time step and stability factor
    tstep.stabilityFactor = 2.0
    dt = 0.8 / 2.0
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Skip adjusting time step
    integrators[1].dt = 2.0
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Skip adjusting time step
    self.assertEqual(dt, tstep.timeStep(integrators))

    # Adjust time step with value bigger than max
    dt = 0.5
    self.assertEqual(dt, tstep.timeStep(integrators))

    return


  def test_currentStep(self):
    """
    Test currentStep().
    """
    tstep = self.tstep
    tstep.maxDt = 10.0
    tstep.dt = tstep.maxDt
    
    self.assertEqual(10.0, tstep.currentStep())

    
    integrators = [Integrator(3.0),
                   Integrator(2.4)]
    dt = 2.4 / 1.2
    tstep.timeStep(integrators)
    self.assertEqual(dt, tstep.currentStep())

    return


# End of file 
