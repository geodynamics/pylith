#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file tests/3dnew/hex8/TestRigidBody.py
##
## @brief Test suite for testing pylith with rigid body motion for
## large deformations in 3-D.

import numpy
from TestHex8 import TestHex8
from rigidbody_soln import AnalyticalSoln
from pylith.utils.VTKDataReader import has_vtk
from pylith.utils.VTKDataReader import VTKDataReader

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class RigidBodyApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="lgdeformrigidbody")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    # Generate spatial databases
    from rigidbody_gendb import GenerateDB
    db = GenerateDB()
    db.run()

    # Run PyLith
    app = RigidBodyApp()
    app.run()
    run_pylith.done = True
  return


class TestRigidBody(TestHex8):
  """
  Test suite for testing pylith with 2-D rigid body motion.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestHex8.setUp(self)
    run_pylith()
    self.outputRoot = "lgdeformrigidbody"
    if has_vtk():
      self.reader = VTKDataReader()
      self.soln = AnalyticalSoln()
    else:
      self.reader = None
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
    pts = numpy.zeros( (ncells, 3), dtype=numpy.float64)
    if name == "total_strain":
      stateVar = self.soln.strain(pts)
    elif name == "stress":
      stateVar = self.soln.stress(pts)
    else:
      raise ValueError("Unknown state variable '%s'." % name)

    return stateVar


# End of file 
