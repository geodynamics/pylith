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

## @file tests/1d/line3/TestAxial.py
##
## @brief Test suite for testing pylith with 1-D axial extension.

import numpy
from TestLine3 import TestLine3
from pylith.utils.VTKDataReader import has_vtk
from pylith.utils.VTKDataReader import VTKDataReader


# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class AxialApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="axialextension")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    app = AxialApp()
    app.run()
    run_pylith.done = True
  return


class TestAxial(TestLine3):
  """
  Test suite for testing pylith with 1-D axial extension.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestLine3.setUp(self)
    run_pylith()
    self.outputRoot = "axialextension"
    if has_vtk():
      self.reader = VTKDataReader()
    else:
      self.reader = None
    return


  def calcDisplacements(self, vertices):
    """
    Calculate displacement field given coordinates of vertices.
    """
    nvertices = self.mesh['nvertices']
    spaceDim = self.mesh['spaceDim']    
    disp = numpy.zeros( (nvertices, spaceDim), dtype=numpy.float64)
    disp[:,0] = -0.2 + 0.1 * vertices[:,0]

    return disp


  def calcStateVar(self, name, vertices, cells):
    """
    Calculate state variable.
    """
    exx = 0.1

    ncells = self.mesh['ncells']
    tensorSize = self.mesh['tensorSize']

    if name == "total_strain":
      stateVar = exx*numpy.ones( (ncells, tensorSize), dtype=numpy.float64)
    
    elif name == "stress":
      lp2m = self.density*self.vp**2
      stateVar = lp2m*exx * numpy.ones( (ncells, tensorSize), 
                                       dtype=numpy.float64)
    else:
      raise ValueError("Unknown state variable '%s'." % name)

    return stateVar


# End of file 
