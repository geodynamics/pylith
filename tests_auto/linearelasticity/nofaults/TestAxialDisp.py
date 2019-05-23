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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/2d/tri3/TestAxialDisp.py
##
## @brief Test suite for testing pylith with 2-D axial extension.

import numpy

from pylith.apps.PyLithApp import PyLithApp
from pylith.tests import run_pylith

from TestTri import TestTri
from axialdisp_soln import AnalyticalSoln
from axialdisp_gendb import GenerateDB

class TestApp(PyLithApp):
  """Local version of PyLithApp.
  """
  def __init__(self):
    PyLithApp.__init__(self, name="axialdisp")
    return


class TestAxialDisp(TestTri):
  """
  Test suite for testing pylith with 2-D axial extension.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestTri.setUp(self)
    run_pylith(TestApp, GenerateDB)
    self.outputRoot = "axialdisp"

    self.soln = AnalyticalSoln()
    return


  def computeVertexField(self, name, vertices):
    """
    Calculate field given coordinates of vertices.
    """
    if name == "displacement":
      return self.soln.displacement(vertices)
    else:
      raise ValueError("Unknown vertex field '{0}'.".format(name))
    return


  def computeCellField(self, name, centroids):
    """
    Calculate field given coordinates of cell centroids.
    """
    if name == "displacement":
      return self.soln.displacement(centroids)
    else:
      raise ValueError("Unknown cell field '{0}'.".format(name))
    return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from .TestAxialDisp import TestAxialDisp as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
