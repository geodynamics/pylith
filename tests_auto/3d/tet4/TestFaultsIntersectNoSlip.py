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

## @file tests/3d/tet4/TestFaultsIntersectNoSlip.py
##
## @brief Test suite for testing pylith with shear slip.

import numpy
from TestTet4 import TestTet4
from sheardisp_soln import AnalyticalSoln

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class FaultsIntersectNoSlipApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="faultsintersectnoslip")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    app = FaultsIntersectNoSlipApp()
    run_pylith.done = True # Put before run() so only called once
    app.run()
  return


class TestFaultsIntersectNoSlip(TestTet4):
  """
  Test suite for testing pylith with shear slip on two faults.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestTet4.setUp(self)
    self.nverticesO = self.mesh['nvertices']

    # Fault x
    self.mesh['nvertices'] += 8
    self.faultMeshX = {'nvertices': 19,
                       'spaceDim': 3,
                       'ncells': 20,
                       'ncorners': 3}

    # Fault y
    self.mesh['nvertices'] += 2
    self.faultMeshY = {'nvertices': 9,
                       'spaceDim': 3,
                       'ncells': 8,
                       'ncorners': 3}
    run_pylith()
    self.outputRoot = "faultsintersectnoslip"

    self.soln = AnalyticalSoln()
    return


  def test_fault_info(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    from pylith.tests.Fault import check_vertex_fields
    fields = ["normal_dir", "final_slip", "slip_time"]

    self.fault = "x"
    filename = "%s-faultx_info.h5" % self.outputRoot
    check_vertex_fields(self, filename, self.faultMeshX, fields)

    self.fault = "y"
    filename = "%s-faulty_info.h5" % self.outputRoot
    check_vertex_fields(self, filename, self.faultMeshY, fields)

    return


  def test_fault_data(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    from pylith.tests.Fault import check_vertex_fields
    fields = ["slip"]

    filename = "%s-faultx.h5" % self.outputRoot
    self.fault = "x"
    check_vertex_fields(self, filename, self.faultMeshX, fields)

    filename = "%s-faulty.h5" % self.outputRoot
    self.fault = "y"
    check_vertex_fields(self, filename, self.faultMeshY, fields)

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

    if self.fault == "x":
      normalDir = (+1.0, 0.0, 0.0)
      finalSlip = 0.0
      faultMesh = self.faultMeshX
    elif self.fault == "y":
      normalDir = (0.0,-1.0, 0.0)
      finalSlip = 0.0
      faultMesh = self.faultMeshY

    slipTime = 0.0
    dim = 3

    nvertices = faultMesh['nvertices']

    if name == "normal_dir":
      field = numpy.zeros( (1, nvertices, dim), dtype=numpy.float64)
      field[0,:,0] = normalDir[0]
      field[0,:,1] = normalDir[1]
      field[0,:,2] = normalDir[2]

    elif name == "final_slip":
      field = numpy.zeros( (1, nvertices, dim), dtype=numpy.float64)
      field[0,:,0] = finalSlip

    elif name == "slip_time":
      field = slipTime*numpy.zeros( (1, nvertices, 1), dtype=numpy.float64)
      
    elif name == "slip":
      field = numpy.zeros( (1, nvertices, dim), dtype=numpy.float64)
      field[0,:,0] = finalSlip

    elif name == "traction_change":
      field = numpy.zeros( (1, nvertices, dim), dtype=numpy.float64)
      field[0,:,0] = 0.0
      
    else:
      raise ValueError("Unknown fault field '%s'." % name)

    return field


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestFaultsIntersectNoSlip import TestFaultsIntersectNoSlip as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
