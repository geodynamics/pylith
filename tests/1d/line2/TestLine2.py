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

## @file tests/1d/line2/TestLine2.py
##
## @brief Generic tests for problems using 1-D bar mesh.

import unittest
import numpy

class TestLine2(unittest.TestCase):
  """
  Generic tests for problems using 1-D bar mesh.
  """

  def setUp(self):
    """
    Setup for tests.
    """
    self.ncells = 4
    self.ncorners = 2
    self.nvertices = 5
    self.spaceDim = 3
    self.tensorSize = 1
    self.vs = 3000.0
    self.vp = 5291.502622129181
    self.density = 2500.0
    return

  def test_elastic_info(self):
    """
    Check elastic info.
    """
    if self.reader is None:
      return

    data = self.reader.read("%s-statevars-elastic_info.vtk" % self.outputRoot)

    # Check cells
    (ncells, ncorners) = data['cells'].shape
    self.assertEqual(self.ncells, ncells)
    self.assertEqual(self.ncorners, ncorners)

    # Check vertices
    (nvertices, spaceDim) = data['vertices'].shape
    self.assertEqual(self.nvertices, nvertices)
    self.assertEqual(self.spaceDim, spaceDim)

    # Check physical properties
    tolerance = 1.0e-5

    # Lame's constant mu (shear modulus)
    muE = self.density*self.vs**2
    diff = numpy.abs(1.0 - data['cell_fields']['mu']/muE)
    okay = diff < tolerance
    if numpy.sum(okay) != ncells:
      print "Expected Lame's constant mu: ",muE
      print "Lame's constant mu: ",data['cell_fields']['mu']
      self.assertEqual(ncells, numpy.sum(okay))    

    # Lame's constant lambda
    lambdaE = self.density*self.vp**2 - 2*muE
    diff = numpy.abs(1.0 - data['cell_fields']['lambda']/lambdaE)
    okay = diff < tolerance
    if numpy.sum(okay) != ncells:
      print "Expected Lame's constant lambda: ",lambdaE
      print "Lame's constant lambda: ",data['cell_fields']['lambda']
      self.assertEqual(ncells, numpy.sum(okay))    

    # Density
    diff = numpy.abs(1.0 - data['cell_fields']['density']/self.density)
    okay = diff < tolerance
    if numpy.sum(okay) != ncells:
      print "Expected density: ",self.density
      print "Density: ",data['cell_fields']['density']
      self.assertEqual(ncells, numpy.sum(okay))    

    return


# End of file
