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

## @file unittests/pytests/problems/TestTimeStepUser.py

## @brief Unit testing of TimeStepUser object.

import unittest
from pylith.problems.TimeStepUser import TimeStepUser

from pyre.units.time import second,year

stepsE = [1.0*year, 2.0*year, 3.0*year]

# ----------------------------------------------------------------------
class TestTimeStepUser(unittest.TestCase):
  """
  Unit testing of TimeStepUser object.
  """

  def setUp(self):
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()
    normalizer.setTimeScale(2.0*second)

    tstep = TimeStepUser()
    tstep._configure()
    tstep.filename = "data/timesteps.txt"
    tstep.preinitialize()
    tstep.initialize(normalizer)
    self.tstep = tstep
    return
  

  def test_constructor(self):
    """
    Test constructor.
    """
    tstep = TimeStepUser()
    tstep._configure()
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    tstep = self.tstep

    for stepE, step in zip(stepsE, tstep.steps):
      valueE = stepE.value / 2.0 # Nondimensionalize
      self.assertEqual(valueE, step)
    return


  def test_numTimeSteps(self):
    """
    Test numTimeSteps().
    """
    tstep = self.tstep

    self.assertEqual(1, tstep.numTimeSteps())

    tstep.totalTime = (12.0*year).value / 2.0 # nondimensionalize
    self.assertEqual(6, tstep.numTimeSteps())

    tstep.loopSteps = True
    tstep.totalTime = (7.0*year).value / 2.0 # nondimensionalize
    self.assertEqual(5, tstep.numTimeSteps())
    return


  def test_timeStep(self):
    """
    Test timeStep().
    """
    tstep = self.tstep

    step1 = (1.0*year).value / 2.0 # nondimensionalize
    step2 = (2.0*year).value / 2.0 # nondimensionalize
    step3 = (3.0*year).value / 2.0 # nondimensionalize

    self.assertEqual(step1, tstep.timeStep(0.5))
    self.assertEqual(step2, tstep.timeStep(0.5))
    self.assertEqual(step3, tstep.timeStep(0.5))
    self.assertEqual(step3, tstep.timeStep(0.5))
    self.assertEqual(step3, tstep.timeStep(0.5))

    tstep.index = 0
    tstep.loopSteps = True
    self.assertEqual(step1, tstep.timeStep(0.5))
    self.assertEqual(step2, tstep.timeStep(0.5))
    self.assertEqual(step3, tstep.timeStep(0.5))
    self.assertEqual(step1, tstep.timeStep(0.5))
    self.assertEqual(step2, tstep.timeStep(0.5))
    return


  def test_currentStep(self):
    """
    Test currentStep().
    """
    tstep = self.tstep

    tstep.timeStep(0.0)
    stepE = 1.0*year
    valueE = stepE.value / 2.0 # nondimensionalize
    self.assertEqual(valueE, tstep.currentStep())
    return


# End of file 
