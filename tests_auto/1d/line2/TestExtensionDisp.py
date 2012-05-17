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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/1d/line2/TestExtensionDisp.py
##
## @brief Test suite for testing pylith with 1-D axial extension.

import numpy

from TestLine2 import TestLine2


# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class AxialApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="extensiondisp")
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


class TestExtensionDisp(TestLine2):
  """
  Test suite for testing pylith with 1-D axial extension.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestLine2.setUp(self)
    run_pylith()
    self.outputRoot = "extensiondisp"
    return


  def calcDisplacements(self, vertices):
    """
    Calculate displacement field given coordinates of vertices.
    """
    nvertices = self.mesh['nvertices']
    spaceDim = self.mesh['spaceDim']    
    disp = numpy.zeros( (1,nvertices,spaceDim), dtype=numpy.float64)
    disp[0,:,0] = -0.2 + 0.1 * vertices[:,0]

    return disp


  def calcStateVar(self, name, vertices, cells):
    """
    Calculate state variable.
    """
    exx = 0.1

    ncells = self.mesh['ncells']
    tensorSize = self.mesh['tensorSize']

    if name == "total_strain":
      stateVar = exx*numpy.ones( (1, ncells, tensorSize), dtype=numpy.float64)
    
    elif name == "stress":
      lp2m = self.density*self.vp**2
      stateVar = lp2m*exx * numpy.ones( (1, ncells, tensorSize), dtype=numpy.float64)
    else:
      raise ValueError("Unknown state variable '%s'." % name)

    return stateVar


# End of file 
