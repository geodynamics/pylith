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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/2d/quad4/TestAxialDisp.py
##
## @brief Test suite for testing pylith with 2-D axial extension.

import numpy
from TestQuad4 import TestQuad4

from axialdisp_soln import AnalyticalSoln

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class AxialApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="axialdisp")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    # Generate spatial databases
    from axialdisp_gendb import GenerateDB
    db = GenerateDB()
    db.run()

    # Run PyLith
    app = AxialApp()
    app.run()
    run_pylith.done = True
  return


class TestAxialDisp(TestQuad4):
  """
  Test suite for testing pylith with 2-D axial extension.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestQuad4.setUp(self)
    run_pylith()
    self.outputRoot = "axialdisp"
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
    elif name == "stress":
      stateVar = self.soln.stress(pts)
    else:
      raise ValueError("Unknown state variable '%s'." % name)

    return stateVar


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestAxialDisp import TestAxialDisp as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
