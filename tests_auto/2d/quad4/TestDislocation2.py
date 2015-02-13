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

## @file tests/2d/quad4/TestDislocation.py
##
## @brief Test suite for testing pylith with slip on two faults.

import numpy
from TestQuad4 import TestQuad4
from dislocation_soln import AnalyticalSoln
from pylith.utils.VTKDataReader import has_vtk
from pylith.utils.VTKDataReader import VTKDataReader
from pylith.tests.Fault import check_vertex_fields

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class LocalApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="dislocation")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    # Run PyLith
    run_pylith.done = True
    app = LocalApp()
    app.run()
  return


class TestDislocation2(TestQuad4):
  """
  Test suite for testing pylith with 2-D axial extension.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestQuad4.setUp(self)
    self.mesh['nvertices'] = 81+18
    self.nverticesO = 81
    self.fault1Mesh = {'nvertices': 9,
                       'spaceDim': 3,
                       'ncells': 8,
                       'ncorners': 2}
    self.fault2Mesh = {'nvertices': 9,
                       'spaceDim': 3,
                       'ncells': 8,
                       'ncorners': 2}

    run_pylith()
    self.outputRoot = "dislocation2"
    if has_vtk():
      self.reader = VTKDataReader()
      self.soln = AnalyticalSoln()
    else:
      self.reader = None

    return


  def test_fault1_info(self):
    """
    Check fault information.
    """
    if self.reader is None:
      return

    filename = "%s-fault1_info.vtk" % self.outputRoot
    fields = ["normal_dir", "final_slip", "slip_time"]
    check_vertex_fields(self, filename, self.fault1Mesh, fields)

    return


  def test_fault1_data(self):
    """
    Check fault information.
    """
    if self.reader is None:
      return

    filename = "%s-fault1_t0000000.vtk" % self.outputRoot
    fields = ["cumulative_slip", "traction_change"]
    check_vertex_fields(self, filename, self.fault1Mesh, fields)

    return


  def test_fault2_info(self):
    """
    Check fault information.
    """
    if self.reader is None:
      return

    filename = "%s-fault2_info.vtk" % self.outputRoot
    fields = ["normal_dir", "final_slip", "slip_time"]
    check_vertex_fields(self, filename, self.fault2Mesh, fields)

    return


  def test_fault2_data(self):
    """
    Check fault information.
    """
    if self.reader is None:
      return

    filename = "%s-fault2_t0000000.vtk" % self.outputRoot
    fields = ["cumulative_slip", "traction_change"]
    check_vertex_fields(self, filename, self.fault2Mesh, fields)

    return


  def calcDisplacements(self, vertices):
    """
    Calculate displacement field given coordinates of vertices.
    """
    return self.soln.displacement(vertices, self.nverticesO)


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

    normalDir = (-1.0, 0.0)
    finalSlip = -2.0
    slipTime = 0.0

    nvertices = self.faultMesh['nvertices']

    if name == "normal_dir":
      field = numpy.zeros( (nvertices, 3), dtype=numpy.float64)
      field[:,0] = normalDir[0]
      field[:,1] = normalDir[1]

    elif name == "final_slip":
      field = numpy.zeros( (nvertices, 3), dtype=numpy.float64)
      field[:,0] = finalSlip
      
    elif name == "slip_time":
      field = slipTime*numpy.zeros( (nvertices, 1), dtype=numpy.float64)
      
    elif name == "cumulative_slip":
      field = numpy.zeros( (nvertices, 3), dtype=numpy.float64)
      field[:,0] = finalSlip

    elif name == "traction_change":
      field = numpy.zeros( (nvertices, 3), dtype=numpy.float64)
      field[:,0] = 0.0
      
    else:
      raise ValueError("Unknown fault field '%s'." % name)

    return field


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestDislocationTwoFaults import TestDislocation2 as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
