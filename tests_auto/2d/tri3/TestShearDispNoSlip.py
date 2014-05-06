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

## @file tests/2d/tri3/TestShearDispNoSlip.py
##
## @brief Test suite for testing pylith with 2-D shear motion with no
## fault slip.

import numpy
from TestTri3 import TestTri3
from sheardisp_soln import AnalyticalSoln

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class ShearApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="sheardispnoslip")
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


class TestShearDispNoSlip(TestTri3):
  """
  Test suite for testing pylith with 2-D shear extension.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestTri3.setUp(self)
    self.nverticesO = self.mesh['nvertices']
    self.mesh['nvertices'] += 1
    self.faultMesh = {'nvertices': 3,
                      'spaceDim': 2,
                      'ncells': 2,
                      'ncorners': 2}
    run_pylith()
    self.outputRoot = "sheardispnoslip"

    self.soln = AnalyticalSoln()
    return


  def test_fault_info(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    filename = "%s-fault_info.h5" % self.outputRoot
    fields = ["normal_dir", "final_slip", "slip_time"]

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

    normalDir = (0.0, 1.0)
    finalSlip = 0.0
    slipTime = 0.0

    nvertices = self.faultMesh['nvertices']

    if name == "normal_dir":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = normalDir[0]
      field[0,:,1] = normalDir[1]

    elif name == "final_slip":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = finalSlip
      
    elif name == "slip_time":
      field = slipTime*numpy.zeros( (1, nvertices, 1), dtype=numpy.float64)
      
    elif name == "slip":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = finalSlip

    elif name == "traction_change":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = 0.0
      
    else:
      raise ValueError("Unknown fault field '%s'." % name)

    # Mask clamped vertices
    maskX = numpy.bitwise_and(vertices[:,0] >= -2.0001e+3, vertices[:,0] <= 1.0)
    maskY = numpy.fabs(vertices[:,1]) <= 0.01e+3
    mask = numpy.bitwise_and(maskX,maskY)
    field[:,~mask,:] = 0

    return field


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestShearDispNoSlip import TestShearDispNoSlip as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
