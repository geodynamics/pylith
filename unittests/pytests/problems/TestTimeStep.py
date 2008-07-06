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

## @file unittests/pytests/problems/TestTimeStep.py

## @brief Unit testing of TimeStep object.

import unittest

from pyre.units.time import second

from pylith.problems.TimeStep import TimeStep

# ----------------------------------------------------------------------
class TestTimeStep(unittest.TestCase):
  """
  Unit testing of TimeStep object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    tstep = TimeStep()
    tstep._configure()
    return
    
  
  def test_preinitialize(self):
    """
    Test preinitialize().
    """
    tstep = TimeStep()
    tstep._configure()
    tstep.preinitialize()
    return


  def test_verifyConfiguration(self):
    """
    Test verifyConfiguration().
    """
    tstep = TimeStep()
    tstep._configure()
    tstep.preinitialize()
    tstep.verifyConfiguration()
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    tstep = TimeStep()
    tstep._configure()
    tstep.preinitialize()
    tstep.initialize()
    return



  def test_numTimeSteps(self):
    """
    Test numTimeSteps().
    """
    tstep = TimeStep()
    tstep._configure()
    try:
      tstep.numTimeSteps()
    except NotImplementedError:
      # do nothing
      x = None
    return
    

  def test_timeStep(self):
    """
    Test timeStep().
    """
    tstep = TimeStep()
    tstep._configure()

    self.assertEqual(0.5*second, tstep.timeStep(0.5*second))
    return


  def test_currentStep(self):
    """
    Test currentStep().
    """
    tstep = TimeStep()
    tstep._configure()

    self.assertEqual(0.0*second, tstep.currentStep())

    tstep.dt = 1.0e-4*second
    self.assertEqual(1.0e-4*second, tstep.currentStep())

    return


# End of file 
