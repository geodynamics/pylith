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

## @file tests/2d/quad4/TestFrictionCompression.py
##
## @brief Test suite for testing pylith with 2-D axial compression with friction.

import numpy
from TestQuad4 import TestQuad4
from friction_compression_soln import AnalyticalSoln
from pylith.utils.VTKDataReader import has_vtk
from pylith.utils.VTKDataReader import VTKDataReader

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class CompressionApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="friction_compression")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    # Run PyLith
    app = CompressionApp()
    app.run()
    run_pylith.done = True
  return


class TestFrictionCompression(TestQuad4):
  """
  Test suite for testing pylith with 2-D axial compression with friction.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestQuad4.setUp(self)
    run_pylith()
    self.outputRoot = "friction_compression"
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
