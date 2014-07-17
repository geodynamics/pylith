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
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/3d/hex8/TestFrictionNoSlip.py
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
    PyLithApp.__init__(self, name="frictionnoslip")
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


class TestFrictionNoSlip(TestHex8):
  """
  Test suite for testing pylith with 2-D shear extension.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestHex8.setUp(self)
    self.nverticesO = self.mesh['nvertices']
    self.mesh['nvertices'] += 10
    self.faultMesh = {'nvertices': 21,
                      'spaceDim': 3,
                      'ncells': 12,
                      'ncorners': 4}

    run_pylith()
    self.outputRoot = "frictionnoslip"

    self.soln = AnalyticalSoln()
    return


  def test_fault_info(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    filename = "%s-fault_info.h5" % self.outputRoot
    fields = ["normal_dir","strike_dir","static_coefficient", "traction_initial"]

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


    if name == "normal_dir":
      field = numpy.zeros( (1, nvertices, 3), dtype=numpy.float64)
      field[0,:,0] = normalDir[0]
      field[0,:,1] = normalDir[1]
      field[0,:,2] = normalDir[2]

    elif name == "strike_dir":
      field = numpy.zeros( (1, nvertices, 3), dtype=numpy.float64)
      field[0,:,0] = strikeDir[0]
      field[0,:,1] = strikeDir[1]
      field[0,:,2] = strikeDir[2]

    elif name == "static_coefficient":
      field = staticCoef*numpy.ones( (1, nvertices, 1), dtype=numpy.float64)
      
    elif name == "traction_initial":
      field = numpy.zeros( (1, nvertices, 3), dtype=numpy.float64)
      field[0,:,0] = initialTraction[0]
      field[0,:,1] = initialTraction[1]
      field[0,:,2] = initialTraction[2]
      
    elif name == "slip":
      field = numpy.zeros( (1, nvertices, 3), dtype=numpy.float64)

    else:
      raise ValueError("Unknown fault field '%s'." % name)

    # Mask clamped vertices
    maskX = numpy.fabs(vertices[:,0]) <= 1.0
    maskY = numpy.fabs(vertices[:,1]) <= 25.001e+3
    maskZ = vertices[:,2] >= -20.001e+3
    mask = numpy.bitwise_and(numpy.bitwise_and(maskX,maskY),maskZ)
    field[:,~mask,:] = 0.0

    return field


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestFrictionNoSlip import TestFrictionNoSlip as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
