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

## @file unittests/pytests/problems/TestTimeStepAdapt.py

## @brief Unit testing of TimeStepAdapt object.

import unittest
from pylith.problems.TimeStepAdapt import TimeStepAdapt

from pyre.units.time import second

# ----------------------------------------------------------------------
class Integrator:

  def __init__(self, dt):
    self.dt = dt


  def stableTimeStep(self, mesh):
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

    self.assertEqual(0.0, tstep.totalTimeN)
    self.assertEqual(0.5, tstep.maxDtN)
    return
  

  def test_numTimeSteps(self):
    """
    Test numTimeSteps().
    """
    tstep = self.tstep

    self.assertEqual(1, tstep.numTimeSteps())

    tstep.totalTimeN = 4.0
    tstep.maxDtN = 2.0
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

    from pylith.topology.Mesh import Mesh
    mesh = Mesh()

    # Set time step
    dt = 0.5 / 2.0
    self.assertEqual(dt, tstep.timeStep(mesh, integrators))

    # Increase stable time step, but time step should not change (skipped)
    integrators[1].dt = 0.8
    self.assertEqual(dt, tstep.timeStep(mesh, integrators))

    # Reduce time step even if though should have skipped
    integrators[1].dt = 0.2
    dt = 0.2 / 2.0
    self.assertEqual(dt, tstep.timeStep(mesh, integrators))

    # Skip adjusting time step
    integrators[1].dt = 0.8
    self.assertEqual(dt, tstep.timeStep(mesh, integrators))

    # Skip adjusting time step
    self.assertEqual(dt, tstep.timeStep(mesh, integrators))

    # Adjust time step and stability factor
    tstep.stabilityFactor = 3.0
    dt = 0.8 / 3.0
    self.assertEqual(dt, tstep.timeStep(mesh, integrators))

    # Skip adjusting time step
    integrators[1].dt = 2.0
    self.assertEqual(dt, tstep.timeStep(mesh, integrators))

    # Skip adjusting time step
    self.assertEqual(dt, tstep.timeStep(mesh, integrators))

    # Adjust time step with value bigger than max
    dt = 0.5
    self.assertEqual(dt, tstep.timeStep(mesh, integrators))

    return


  def test_currentStep(self):
    """
    Test currentStep().
    """
    tstep = self.tstep
    tstep.maxDtN = 10.0
    tstep.dtN = tstep.maxDtN
    
    self.assertEqual(10.0, tstep.currentStep())

    
    integrators = [Integrator(3.0),
                   Integrator(2.4)]

    from pylith.topology.Mesh import Mesh
    mesh = Mesh()

    dt = 2.4 / 2.0
    tstep.timeStep(mesh, integrators)
    self.assertEqual(dt, tstep.currentStep())

    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.problems.TimeStepAdapt import time_step
    ts = time_step()
    return


# End of file 
