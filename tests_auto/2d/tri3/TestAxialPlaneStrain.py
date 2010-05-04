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

## @file tests/2d/tri3/TestAxialPlaneStrain.py
##
## @brief Test suite for testing pylith with axial compression in
## y-direction for 2-D box.

import unittest
import numpy
from pylith.utils.VTKDataReader import has_vtk
from pylith.utils.VTKDataReader import VTKDataReader

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class AxialPlaneStrainApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="axialplanestrain")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    app = AxialPlaneStrainApp()
    app.run()
    run_pylith.done = True
  return


class TestAxialPlaneStrain(unittest.TestCase):
  """
  Test suite for testing pylith with axial extension in y-direction
  for 2-D box.
  """

  def setUp(self):
    """
    Setup for test.
    """
    run_pylith()
    if has_vtk():
      self.reader = VTKDataReader()
    else:
      self.reader = None
    return

  
  def test_elastic_info(self):
    """
    Check elastic info.
    """
    if self.reader is None:
      return

    data = self.reader.read("axialplanestrain-statevars-elastic_info.vtk")

    # Check cells
    ncellsE = 84
    ncornersE = 3
    (ncells, ncorners) = data['cells'].shape
    self.assertEqual(ncellsE, ncells)
    self.assertEqual(ncornersE, ncorners)

    # Check vertices
    nverticesE = 54
    spaceDimE = 3
    (nvertices, spaceDim) = data['vertices'].shape
    self.assertEqual(nverticesE, nvertices)
    self.assertEqual(spaceDimE, spaceDim)

    # Check physical properties
    tolerance = 1.0e-5
    vsE = 3000.0
    vpE = 5291.502622129181
    densityE = 2500.0

    # Lame's constant mu (shear modulus)
    muE = densityE*vsE**2
    diff = numpy.abs(1.0 - data['cell_fields']['mu']/muE)
    okay = diff < tolerance
    if numpy.sum(okay) != ncells:
      print "Lame's constant mu: ",data['cell_fields']['mu']
      self.assertEqual(ncells, numpy.sum(okay))    

    # Lame's constant lambda
    lambdaE = densityE*vpE**2 - 2*muE
    diff = numpy.abs(1.0 - data['cell_fields']['lambda']/lambdaE)
    okay = diff < tolerance
    if numpy.sum(okay) != ncells:
      print "Lame's constant lambda: ",data['cell_fields']['lambda']
      self.assertEqual(ncells, numpy.sum(okay))    

    # Density
    diff = numpy.abs(1.0 - data['cell_fields']['density']/densityE)
    okay = diff < tolerance
    if numpy.sum(okay) != ncells:
      print "Density: ",data['cell_fields']['density']
      self.assertEqual(ncells, numpy.sum(okay))    
    return


  def test_soln(self):
    """
    Check solution (displacement) field.
    """
    if self.reader is None:
      return

    data = self.reader.read("axialplanestrain_t0000000.vtk")

    # Check cells
    ncellsE = 84
    ncornersE = 3
    (ncells, ncorners) = data['cells'].shape
    self.assertEqual(ncellsE, ncells)
    self.assertEqual(ncornersE, ncorners)

    # Check vertices
    nverticesE = 54
    spaceDimE = 3
    vertices = data['vertices']
    (nvertices, spaceDim) = vertices.shape
    self.assertEqual(nverticesE, nvertices)
    self.assertEqual(spaceDimE, spaceDim)

    # Check displacement solution
    tolerance = 1.0e-5
    dispE = numpy.zeros( (nvertices, spaceDim), dtype=numpy.float64)
    dispE[:,1] = -0.004 * vertices[:,1]

    disp = data['vertex_fields']['displacements']

    # Check x displacements
    diff = numpy.abs(disp[:,0] - dispE[:,0])
    okay = diff < tolerance
    if numpy.sum(okay) != nvertices:
      print "Displacement field: ",disp
      self.assertEqual(nvertices, numpy.sum(okay))    
    
    # Check y displacements
    mask = dispE[:,1] > 0.0
    diff = mask * numpy.abs(1.0 - disp[:,1] / dispE[:,1]) + \
        ~mask * numpy.abs(disp[:,1] - dispE[:,1])
    okay = diff < tolerance
    if numpy.sum(okay) != nvertices:
      print "Displacement field: ",disp
      self.assertEqual(nvertices, numpy.sum(okay))    

    # Check z displacements
    diff = numpy.abs(disp[:,2] - dispE[:,2])
    okay = diff < tolerance
    if numpy.sum(okay) != nvertices:
      print "Displacement field: ",disp
      self.assertEqual(nvertices, numpy.sum(okay))    
    
    return


  def test_elastic_statevars(self):
    """
    Check elastic state variables.
    """
    return


# End of file 
