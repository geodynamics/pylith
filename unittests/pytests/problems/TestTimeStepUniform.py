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

from spatialdata.units.Nondimensional import Nondimensional
from pyre.units.time import second

# ----------------------------------------------------------------------
class TestTimeStepUniform(unittest.TestCase):
  """
  Unit testing of TimeStepUniform object.
  """

  def setUp(self):
    """
    Setup time step object.
    """
    normalizer = Nondimensional()
    normalizer._configure()

    tstep = TimeStepUniform()
    tstep._configure()
    tstep.preinitialize()
    tstep.verifyConfiguration()
    tstep.initialize(normalizer)
    self.tstep = tstep
    return


  def test_numTimeSteps(self):
    """
    Test numTimeSteps().
    """
    tstep = self.tstep

    self.assertEqual(1, tstep.numTimeSteps())

    tstep.totalTimeN = 4.0
    tstep.dtN = 2.0
    self.assertEqual(3, tstep.numTimeSteps())

    return


  def test_timeStep(self):
    """
    Test timeStep().
    """
    tstep = self.tstep

    integrators = None
    mesh = None

    self.assertEqual(1.0, tstep.timeStep(mesh, integrators))

    tstep.dtN = 1.0e-4
    self.assertEqual(1.0e-4, tstep.timeStep(mesh, integrators))

    return


  def test_currentStep(self):
    """
    Test currentStep().
    """
    tstep = self.tstep

    self.assertEqual(1.0, tstep.currentStep())

    tstep.dtN = 1.0e-4
    self.assertEqual(1.0e-4, tstep.currentStep())

    return


# End of file 
