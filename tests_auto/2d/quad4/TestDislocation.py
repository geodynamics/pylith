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
## @brief Test suite for testing pylith with fault slip.

import numpy
from TestQuad4 import TestQuad4
from dislocation_soln import AnalyticalSoln

from pylith.tests.Fault import check_vertex_fields

# ----------------------------------------------------------------------
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


# ----------------------------------------------------------------------
class TestDislocation(TestQuad4):
  """
  Test suite for fault with prescribed slip.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestQuad4.setUp(self)
    self.mesh['nvertices'] = 81+9
    self.nverticesO = 81
    self.faultMesh = {'nvertices': 9,
                      'spaceDim': 2,
                      'ncells': 8,
                      'ncorners': 2}

    run_pylith()
    self.outputRoot = "dislocation"
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
    check_vertex_fields(self, filename, self.faultMesh, fields)

    return


  def test_fault_data(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    filename = "%s-fault.h5" % self.outputRoot
    fields = ["slip", "traction_change"]
    check_vertex_fields(self, filename, self.faultMesh, fields)

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

    return field


# ----------------------------------------------------------------------
# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class LocalApp2(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="dislocation_np2")
    return


# Helper function to run PyLith
def run_pylith2():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith2):
    # Run PyLith
    run_pylith2.done = True
    app = LocalApp2()
    app.run()
  return


# ----------------------------------------------------------------------
class TestDislocation2(TestDislocation):
  """
  Test suite for fault with prescribed slip w/2 procs.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestQuad4.setUp(self)
    self.mesh['nvertices'] = 81+9
    self.nverticesO = 81
    self.faultMesh = {'nvertices': 9,
                      'spaceDim': 2,
                      'ncells': 8,
                      'ncorners': 2}

    run_pylith2()
    self.outputRoot = "dislocation_np2"
    self.soln = AnalyticalSoln()

    return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestDislocation import TestDislocation as Tester
  from TestDislocation import TestDislocation2 as Tester2

  suite = unittest.TestSuite()

  suite.addTest(unittest.makeSuite(Tester))
  suite.addTest(unittest.makeSuite(Tester2))

  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
