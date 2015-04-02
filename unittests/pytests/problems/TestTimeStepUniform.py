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

## @file unittests/pytests/problems/TestTimeStepUniform.py

## @brief Unit testing of TimeStepUniform object.

import unittest
from pylith.problems.TimeStepUniform import TimeStepUniform

from spatialdata.units.Nondimensional import Nondimensional
from pyre.units.time import second

# ----------------------------------------------------------------------
class Integrator:

  def __init__(self, dt):
    self.dt = dt
    return


  def stableTimeStep(self, mesh):
    return self.dt


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

    integrators = [Integrator(4.0),
                   Integrator(8.0)]

    from pylith.topology.Mesh import Mesh
    mesh = Mesh()

    self.assertEqual(1.0, tstep.timeStep(mesh, integrators))

    tstep.dtN = 0.5
    self.assertEqual(0.5, tstep.timeStep(mesh, integrators))

    caught = False
    try:
      tstep.dtN = 10.0
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

    self.assertEqual(1.0, tstep.currentStep())

    tstep.dtN = 1.0e-4
    self.assertEqual(1.0e-4, tstep.currentStep())

    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.problems.TimeStepUniform import time_step
    ts = time_step()
    return


# End of file 
