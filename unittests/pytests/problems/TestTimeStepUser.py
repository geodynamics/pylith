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
    tstep = TimeStepUser()
    tstep._configure()

    tstep.filename = "data/timesteps.txt"
    tstep.preinitialize()
    tstep.initialize()

    for stepE, step in zip(stepsE, tstep.steps):
      self.assertEqual(stepE, step)
    return


  def test_numTimeSteps(self):
    """
    Test numTimeSteps().
    """
    tstep = TimeStepUser()
    tstep._configure()

    tstep.filename = "data/timesteps.txt"
    tstep.preinitialize()
    tstep.initialize()

    self.assertEqual(1, tstep.numTimeSteps())

    tstep.totalTime = 7.0*year
    self.assertEqual(5, tstep.numTimeSteps())
    return


  def test_timeStep(self):
    """
    Test timeStep().
    """
    tstep = TimeStepUser()
    tstep._configure()

    tstep.filename = "data/timesteps.txt"
    tstep.preinitialize()
    tstep.initialize()

    self.assertEqual(1.0*year, tstep.timeStep(0.5*year))
    self.assertEqual(2.0*year, tstep.timeStep(0.5*year))
    self.assertEqual(3.0*year, tstep.timeStep(0.5*year))
    self.assertEqual(1.0*year, tstep.timeStep(0.5*year))
    self.assertEqual(2.0*year, tstep.timeStep(0.5*year))
    return


  def test_currentStep(self):
    """
    Test currentStep().
    """
    tstep = TimeStepUser()
    tstep._configure()

    tstep.filename = "data/timesteps.txt"
    tstep.preinitialize()
    tstep.initialize()

    tstep.timeStep(0.0*second)
    self.assertEqual(1.0*year, tstep.currentStep())
    return


# End of file 
