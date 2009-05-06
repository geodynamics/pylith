#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
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
    self.nvertices = 6
    self.nverticesO = 5

    run_pylith()
    self.outputRoot = "dislocation"
    if has_vtk():
      self.reader = VTKDataReader()
    else:
      self.reader = None
    return


  def test_soln(self):
    """
    Check solution (displacement) field.
    """
    if self.reader is None:
      return

    data = self.reader.read("%s_t0000000.vtk" % self.outputRoot)

    # Check cells
    (ncells, ncorners) = data['cells'].shape
    self.assertEqual(self.ncells, ncells)
    self.assertEqual(self.ncorners, ncorners)

    # Check vertices
    vertices = data['vertices']
    (nvertices, spaceDim) = vertices.shape
    self.assertEqual(self.nvertices, nvertices)
    self.assertEqual(self.spaceDim, spaceDim)

    # Check displacement solution
    tolerance = 1.0e-5
    dispE = numpy.zeros( (nvertices, spaceDim), dtype=numpy.float64)
    maskP = vertices[:,0] >= 2.0
    maskP[self.nverticesO:self.nvertices] = False
    maskN = numpy.bitwise_and(vertices[:,0] <= 2.0, ~maskP)
    dispE[:,0] = \
        maskN*(-0.20 - 0.025*vertices[:,0]) + \
        maskP*(+0.30 - 0.025*vertices[:,0])

    disp = data['vertex_fields']['displacement']

    # Check x displacements
    diff = numpy.abs(disp[:,0] - dispE[:,0])
    okay = diff < tolerance
    if numpy.sum(okay) != nvertices:
      print "Displacement field expected: ",dispE
      print "Displacement field: ",disp
    self.assertEqual(nvertices, numpy.sum(okay))    
    
    # Check y displacements
    diff = numpy.abs(disp[:,1] - dispE[:,1])
    okay = diff < tolerance
    if numpy.sum(okay) != nvertices:
      print "Displacement field expected: ",dispE
      print "Displacement field: ",disp
    self.assertEqual(nvertices, numpy.sum(okay))    

    # Check z displacements
    diff = numpy.abs(disp[:,2] - dispE[:,2])
    okay = diff < tolerance
    if numpy.sum(okay) != nvertices:
      print "Displacement field expected: ",dispE
      print "Displacement field: ",disp
    self.assertEqual(nvertices, numpy.sum(okay))    
    
    return


  def test_elastic_statevars(self):
    """
    Check elastic state variables.
    """
    if self.reader is None:
      return

    data = self.reader.read("%s-statevars-elastic_t0000000.vtk" % \
                              self.outputRoot)

    # Check cells
    (ncells, ncorners) = data['cells'].shape
    self.assertEqual(self.ncells, ncells)
    self.assertEqual(self.ncorners, ncorners)

    # Check vertices
    vertices = data['vertices']
    (nvertices, spaceDim) = vertices.shape
    self.assertEqual(self.nvertices, nvertices)
    self.assertEqual(self.spaceDim, spaceDim)

    # Check strains
    tolerance = 1.0e-5
    exx = -0.025
    strainE = exx*numpy.ones( (ncells, self.tensorSize), dtype=numpy.float64)
    strain = data['cell_fields']['total_strain']
    diff = numpy.abs(strain[:]-strainE[:,0])
    okay = diff < tolerance
    if numpy.sum(okay) != ncells:
      print "Strain field expected: ",strainE
      print "Strain field: ",strain
    self.assertEqual(ncells, numpy.sum(okay))    

    # Check stresses
    lp2m = self.density*self.vp**2
    stressE = lp2m*exx * numpy.ones( (ncells, self.tensorSize), 
                                     dtype=numpy.float64)
    stress = data['cell_fields']['stress']
    diff = numpy.abs(stress[:]-stressE[:,0])
    okay = diff < tolerance
    if numpy.sum(okay) != ncells:
      print "Stress field expected: ",stressE
      print "Stress field: ",stress
    self.assertEqual(ncells, numpy.sum(okay))    

    return


# End of file 
