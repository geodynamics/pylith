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

## @file tests/2d/quad4/TestLgDeformTraction.py
##
## @brief Test suite for testing pylith with uniform tractions for
## large deformations in 2-D.

import numpy
from TestQuad4 import TestQuad4

from lgdeformtraction_soln import AnalyticalSoln

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class TractionApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="lgdeformtraction")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    # Run PyLith
    app = TractionApp()
    app.run()
    run_pylith.done = True
  return


class TestTraction(TestQuad4):
  """
  Test suite for testing pylith with 2-D axial tractions.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestQuad4.setUp(self)
    run_pylith()
    self.outputRoot = "lgdeformtraction"
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
  from TestLgDeformTraction import TestLgDeformTraction as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
