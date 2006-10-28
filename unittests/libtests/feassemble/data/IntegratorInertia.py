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

## @file unittests/libtests/feassemble/data/IntegratorInertia.py

## @brief Python application for generating C++ data files for testing
## C++ IntegratorInertia object.

from IntegratorApp import IntegratorApp

import numpy

# ----------------------------------------------------------------------

# IntegratorInertia class
class IntegratorInertia(IntegratorApp):
  """
  Python application for generating C++ data files for testing C++
  IntegratorInertia object.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="itnegratorinertia"):
    """
    Constructor.
    """
    IntegratorApp.__init__(self, name)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _calculateMatrix(self):
    """
    Calculate matrix associated with integration.
    """
    import feutils
    
    self.valsMatrix = numpy.zeros( (self.spaceDim*self.numVertices,
                                    self.spaceDim*self.numVertices),
                                   dtype=numpy.float64)

    n = numpy.zeros( (self.spaceDim, self.spaceDim*self.numCorners),
                     dtype=numpy.float64)
    
    for cell in self.cells:
      cellMatrix = numpy.zeros( (self.spaceDim*self.numCorners,
                                 self.spaceDim*self.numCorners),
                                dtype=numpy.float64)

      vertices = self.vertices[cell, :]
      (jacobian, jacobianInv, jacobianDet) = \
                 feutils.calculateJacobian(self.quadrature, vertices)
      density = 1.0
      for iQuad in xrange(self.numQuadPts):
        n *= 0.0
        for iBasis in xrange(self.numCorners):
          for iDim in xrange(self.spaceDim):
            n[iDim, iBasis*self.spaceDim+iDim] = self.basis[iQuad, iBasis]

        wt = density * self.quadWts[iQuad] * jacobianDet[iQuad]
        cellMatrix[:] += wt * numpy.dot(n.transpose(), n)
      feutils.assembleMat(self.valsMatrix, cellMatrix, cell)
    return


# End of file 
