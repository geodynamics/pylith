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

## @file tests/1d/line2/TestDislocation.py
##
## @brief Test suite for testing pylith with 1-D axial extension.

import numpy
from TestLine2 import TestLine2
from pylith.utils.VTKDataReader import has_vtk
from pylith.utils.VTKDataReader import VTKDataReader
from pylith.tests.Fault import check_vertex_fields

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class DislocationApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="dislocation")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    app = DislocationApp()
    app.run()
    run_pylith.done = True
  return


class TestDislocation(TestLine2):
  """
  Test suite for testing pylith with 1-D axial extension.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestLine2.setUp(self)
    self.mesh['nvertices'] = 6
    self.nverticesO = 5

    self.faultMesh = {'nvertices': 1,
                      'spaceDim': 3,
                      'ncells': 1,
                      'ncorners': 1}

    run_pylith()
    self.outputRoot = "dislocation"
    if has_vtk():
      self.reader = VTKDataReader()
    else:
      self.reader = None
    return


  def test_fault_info(self):
    """
    Check fault information.
    """
    if self.reader is None:
      return

    filename = "%s-fault_info.vtk" % self.outputRoot
    fields = ["normal_dir", "final_slip", "slip_time"]
    check_vertex_fields(self, filename, self.faultMesh, fields)

    return


  def test_fault_data(self):
    """
    Check fault information.
    """
    if self.reader is None:
      return

    filename = "%s-fault_t0000000.vtk" % self.outputRoot
    fields = ["slip", "traction_change"]
    check_vertex_fields(self, filename, self.faultMesh, fields)

    return


  def calcDisplacements(self, vertices):
    """
    Calculate displacement field given coordinates of vertices.
    """
    nvertices = self.mesh['nvertices']
    spaceDim = self.mesh['spaceDim']    
    nverticesO = self.nverticesO

    disp = numpy.zeros( (nvertices, spaceDim), dtype=numpy.float64)
    maskP = vertices[:,0] >= 2.0
    maskP[nverticesO:nvertices] = False
    maskN = numpy.bitwise_and(vertices[:,0] <= 2.0, ~maskP)
    disp[:,0] = \
        maskN*(-0.20 - 0.025*vertices[:,0]) + \
        maskP*(+0.30 - 0.025*vertices[:,0])

    return disp


  def calcStateVar(self, name, vertices, cells):
    """
    Calculate state variable.
    """
    exx = -0.025

    ncells = self.mesh['ncells']
    tensorSize = self.mesh['tensorSize']

    if name == "total_strain":
      stateVar = exx*numpy.ones( (ncells, tensorSize), dtype=numpy.float64)
    
    elif name == "stress":
      lp2m = self.density*self.vp**2
      stateVar = lp2m*exx * numpy.ones( (ncells, tensorSize), 
                                       dtype=numpy.float64)
    else:
      raise ValueError("Unknown state variable '%s'." % name)

    return stateVar


  def calcFaultField(self, name, vertices):
    """
    Calculate fault info.
    """

    normalDir = 1.0
    finalSlip = -0.5
    slipTime = 0.0

    slip = 1.0
    exx = -0.025
    lp2m = self.density*self.vp**2
    traction = -exx*lp2m

    nvertices = self.faultMesh['nvertices']

    if name == "normal_dir":
      field = numpy.zeros( (nvertices, 3), dtype=numpy.float64)
      field[:,0] = normalDir

    elif name == "final_slip":
      field = numpy.zeros( (nvertices, 3), dtype=numpy.float64)
      field[:,0] = finalSlip
      
    elif name == "slip_time":
      field = slipTime*numpy.ones( (nvertices, 1), dtype=numpy.float64)
      
    elif name == "slip":
      field = numpy.zeros( (nvertices, 3), dtype=numpy.float64)
      field[:,0] = finalSlip

    elif name == "traction_change":
      field = numpy.zeros( (nvertices, 3), dtype=numpy.float64)
      field[:,0] = traction
      
    else:
      raise ValueError("Unknown fault field '%s'." % name)

    return field


# End of file 
