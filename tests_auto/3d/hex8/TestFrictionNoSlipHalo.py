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

## @file tests/3d/hex8/TestFrictionNoSlipHalo.py
##
## @brief Test suite for testing pylith with 3-D shear motion with no
## fault slip.

import numpy
from TestHex8 import TestHex8
from sheardisp_soln import AnalyticalSoln

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class ShearApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="frictionnoslip_halo")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    # Generate spatial databases
    from sheardisp_gendb import GenerateDB
    db = GenerateDB()
    db.run()

    # Run PyLith
    app = ShearApp()
    run_pylith.done = True # Put before run() so only called once
    app.run()
  return


class TestFrictionNoSlipHalo(TestHex8):
  """
  Test suite for testing pylith with 2-D shear extension.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestHex8.setUp(self)
    self.nverticesO = self.mesh['nvertices']
    self.mesh['nvertices'] += 21
    self.faultMesh = {'nvertices': 34,
                      'spaceDim': 3,
                      'ncells': 23,
                      'ncorners': 4}

    run_pylith()
    self.outputRoot = "frictionnoslip_halo"

    self.soln = AnalyticalSoln()
    return


  def test_fault_info(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    filename = "%s-fault_info.h5" % self.outputRoot
    fields = ["static_coefficient"]

    from pylith.tests.Fault import check_vertex_fields
    check_vertex_fields(self, filename, self.faultMesh, fields)

    return


  def test_fault_data(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    filename = "%s-fault.h5" % self.outputRoot
    fields = ["slip"]

    from pylith.tests.Fault import check_vertex_fields
    check_vertex_fields(self, filename, self.faultMesh, fields)

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


  def calcFaultField(self, name, vertices):
    """
    Calculate fault info.
    """

    normalDir = (-1.0, 0.0, 0.0)
    strikeDir = (0.0, -1.0, 0.0)

    staticCoef  = 0.6
    initialTraction = (0.0,0.0,-100.0e+6)

    nvertices = self.faultMesh['nvertices']


    if name == "static_coefficient":
      field = staticCoef*numpy.ones( (1, nvertices, 1), dtype=numpy.float64)
      
    elif name == "slip":
      field = numpy.zeros( (1, nvertices, 3), dtype=numpy.float64)

    else:
      raise ValueError("Unknown fault field '%s'." % name)

    return field


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestFrictionNoSlipHalo import TestFrictionNoSlipHalo as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
