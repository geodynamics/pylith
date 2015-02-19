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

## @file unittests/pytests/problems/TestTimeStepUser.py

## @brief Unit testing of TimeStepUser object.

import unittest
from pylith.problems.TimeStepUser import TimeStepUser

from pyre.units.time import second,year

stepsE = [2*1.0, 2*2.0, 2*3.0]

# ----------------------------------------------------------------------
class Integrator:

  def __init__(self, dt):
    self.dt = dt


  def stableTimeStep(self, mesh):
    return self.dt


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

    integrators = [Integrator(40.0),
                   Integrator(80.0)]

    from pylith.topology.Mesh import Mesh
    mesh = Mesh()

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

    integrators = [Integrator(0.01),
                   Integrator(8.0)]
    caught = False
    try:
      tstep.timeStep(mesh, integrators)
    except RuntimeError:
      caught = True
    self.failUnless(caught)

    return


  def test_currentStep(self):
    """
    Test currentStep().
    """
    tstep = self.tstep

    integrators = [Integrator(4.0),
                   Integrator(8.0)]

    from pylith.topology.Mesh import Mesh
    from pylith.mpi.Communicator import petsc_comm_world
    mesh = Mesh()
    #mesh.setComm(petsc_comm_world())

    tstep.timeStep(mesh, integrators)
    stepE = 1.0 / 0.5 # Nondimensionalize
    self.assertEqual(stepE, tstep.currentStep())
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.problems.TimeStepUser import time_step
    ts = time_step()
    return


# End of file 
