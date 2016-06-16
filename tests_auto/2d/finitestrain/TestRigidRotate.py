#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/2d/quad4/TestRigidRotate.py
##
## @brief Test suite for testing pylith with 2-D rigid body rotation.

import numpy
from TestQuad4 import TestQuad4

from rigidrotate_soln import AnalyticalSoln,p_mu

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class RotateApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="rigidrotate")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    # Generate spatial databases
    from rigidrotate_gendb import GenerateDB
    db = GenerateDB()
    db.run()

    # Run PyLith
    app = RotateApp()
    app.run()
    run_pylith.done = True
  return


class TestRigidRotate(TestQuad4):
  """
  Test suite for testing pylith with 2-D rigid body rotation.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestQuad4.setUp(self)
    run_pylith()
    self.outputRoot = "rigidrotate"
    self.soln = AnalyticalSoln()
    return


  def calcDisplacements(self, vertices):
    """
    Calculate displacement field given coordinates of vertices.
    """
    return self.soln.displacement(vertices)


  def calcStateVar(self, name, vertices, cells):
    """
    Calculate state variable.
    """
    ncells = self.mesh['ncells']
    pts = numpy.zeros( (ncells, 2), dtype=numpy.float64)
    if name == "total_strain":
      stateVar = self.soln.strain(pts)
    elif name == "stress" or name == "cauchy_stress":
      stateVar = self.soln.stress(pts)
    else:
      raise ValueError("Unknown state variable '%s'." % name)

    return stateVar


  def getValueScale(self, name):
    """
    Get scale for value.
    """
    if name == "total_strain":
      scale = 1.0
    elif name == "stress":
      scale = p_mu
    elif name == "cauchy_stress":
      scale = p_mu
    else:
      raise ValueError("Unknown variable '%s'." % name)

    return scale


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestRigidRotate import TestRigidRotate as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
