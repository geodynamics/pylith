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

stepsE = [2*1.0, 2*2.0, 2*3.0]

# ----------------------------------------------------------------------
class TestTimeStepUser(unittest.TestCase):
  """
  Unit testing of TimeStepUser object.
  """

  def setUp(self):
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()
    normalizer.setTimeScale(0.5*year)

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
      self.assertEqual(stepE, step)
    return


  def test_numTimeSteps(self):
    """
    Test numTimeSteps().
    """
    tstep = self.tstep

    self.assertEqual(1, tstep.numTimeSteps())

    tstep.totalTimeN = 12.0 / 0.5 # nondimensionalize
    self.assertEqual(6, tstep.numTimeSteps())

    tstep.loopSteps = True
    tstep.totalTimeN = 7.0 / 0.5 # nondimensionalize
    self.assertEqual(5, tstep.numTimeSteps())
    return


  def test_timeStep(self):
    """
    Test timeStep().
    """
    tstep = self.tstep

    step1 = 1.0 / 0.5 # nondimensionalize
    step2 = 2.0 / 0.5 # nondimensionalize
    step3 = 3.0 / 0.5 # nondimensionalize

    integrators = None
    mesh = None

    self.assertEqual(step1, tstep.timeStep(mesh, integrators))
    self.assertEqual(step2, tstep.timeStep(mesh, integrators))
    self.assertEqual(step3, tstep.timeStep(mesh, integrators))
    self.assertEqual(step3, tstep.timeStep(mesh, integrators))
    self.assertEqual(step3, tstep.timeStep(mesh, integrators))

    tstep.index = 0
    tstep.loopSteps = True
    self.assertEqual(step1, tstep.timeStep(mesh, integrators))
    self.assertEqual(step2, tstep.timeStep(mesh, integrators))
    self.assertEqual(step3, tstep.timeStep(mesh, integrators))
    self.assertEqual(step1, tstep.timeStep(mesh, integrators))
    self.assertEqual(step2, tstep.timeStep(mesh, integrators))
    return


  def test_currentStep(self):
    """
    Test currentStep().
    """
    tstep = self.tstep

    integrators = None
    mesh = None

    tstep.timeStep(mesh, integrators)
    stepE = 1.0 / 0.5 # Nondimensionalize
    self.assertEqual(stepE, tstep.currentStep())
    return


# End of file 
