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

## @file tests/2d/quad4/TestFrictionCompression.py
##
## @brief Test suite for testing pylith with 2-D axial compression with friction.

import numpy
from TestQuad4 import TestQuad4
from friction_compression_soln import AnalyticalSoln

from pylith.tests.Fault import check_vertex_fields

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class CompressionApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="friction_compression")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    # Run PyLith
    app = CompressionApp()
    app.run()
    run_pylith.done = True
  return


class TestFrictionCompression(TestQuad4):
  """
  Test suite for testing pylith with 2-D axial compression with friction.
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
    self.outputRoot = "friction_compression"
    self.soln = AnalyticalSoln()
    return


  def test_fault_info(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    filename = "%s-fault_info.h5" % self.outputRoot
    fields = ["strike_dir", "normal_dir", "traction_initial","friction_coefficient","cohesion"]
    check_vertex_fields(self, filename, self.faultMesh, fields)

    return


  def test_fault_data(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    filename = "%s-fault.h5" % self.outputRoot
    fields = ["slip", "traction"]
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

    strikeDir = (0.0, -1.0)
    normalDir = (-1.0, 0.0)
    initialTraction = (0.0, -1.0e+6)
    frictionCoefficient = 0.6

    nvertices = self.faultMesh['nvertices']

    if name == "strike_dir":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = strikeDir[0]
      field[0,:,1] = strikeDir[1]

    elif name == "normal_dir":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = normalDir[0]
      field[0,:,1] = normalDir[1]

    elif name == "traction_initial":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = initialTraction[0]
      field[0,:,1] = initialTraction[1]

    elif name == "friction_coefficient":
      field = numpy.zeros( (1, nvertices, 1), dtype=numpy.float64)
      field[0,:] = frictionCoefficient

    elif name == "cohesion":
      field = numpy.zeros( (1, nvertices, 1), dtype=numpy.float64)

    elif name == "slip":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)

    elif name == "traction":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = 0.0
      field[0,:,1] = -2.2e+6
      
    else:
      raise ValueError("Unknown fault field '%s'." % name)

    return field


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestFrictionCompression import TestFrictionCompression as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
