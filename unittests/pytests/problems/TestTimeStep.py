#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
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

  def setUp(self):
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()
    normalizer.setTimeScale(2.0*second)

    tstep = TimeStep()
    tstep._configure()
    tstep.preinitialize()
    tstep.initialize(normalizer)
    self.tstep = tstep
    return
  

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
    tstep = self.tstep
    return



  def test_numTimeSteps(self):
    """
    Test numTimeSteps().
    """
    tstep = self.tstep
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
    tstep = self.tstep
    tstep._configure()

    integrators = None
    mesh = None
    self.assertEqual(0.0, tstep.timeStep(mesh, integrators))
    return


  def test_currentStep(self):
    """
    Test currentStep().
    """
    tstep = self.tstep

    self.assertEqual(0.0, tstep.currentStep())

    tstep.dtN = 1.0e-4
    self.assertEqual(1.0e-4, tstep.currentStep())

    return


# End of file 
