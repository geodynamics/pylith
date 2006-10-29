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
  
  def __init__(self, name="integratorinertia"):
    """
    Constructor.
    """
    IntegratorApp.__init__(self, name)

    self.valsLumped = None
    return
  

  def _initData(self):
    IntegratorApp._initData(self)
    self.data.addArray(vtype="double", name="_valsLumped",
                       values=self.valsLumped,
                       format="%16.8e", ncols=self.spaceDim)
    return

    
  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _calculateMatrix(self):
    """
    Calculate matrix associated with integration.
    """
    import feutils
    
    self.valsMatrix = numpy.zeros( (self.fiberDim*self.numVertices,
                                    self.fiberDim*self.numVertices),
                                   dtype=numpy.float64)
    self.valsLumped = numpy.zeros( (self.fiberDim*self.numVertices,),
                                   dtype=numpy.float64)

    n = numpy.zeros( (self.fiberDim, self.fiberDim*self.numCorners),
                     dtype=numpy.float64)
    
    for cell in self.cells:
      cellMatrix = numpy.zeros( (self.fiberDim*self.numCorners,
                                 self.fiberDim*self.numCorners),
                                dtype=numpy.float64)

      vertices = self.vertices[cell, :]
      (jacobian, jacobianInv, jacobianDet) = \
                 feutils.calculateJacobian(self.quadrature, vertices)
      density = 1.0
      for iQuad in xrange(self.numQuadPts):
        n *= 0.0
        for iBasis in xrange(self.numCorners):
          for iDim in xrange(self.fiberDim):
            n[iDim, iBasis*self.fiberDim+iDim] = self.basis[iQuad, iBasis]

        wt = density * self.quadWts[iQuad] * jacobianDet[iQuad]
        cellMatrix[:] += wt * numpy.dot(n.transpose(), n)
      feutils.assembleMat(self.valsMatrix, cellMatrix, cell, self.fiberDim)

      cellVec = self._calculateMatrixLumped(cellMatrix)
      feutils.assembleVec(self.valsLumped, cellVec, cell, self.fiberDim)
    return


  def _calculateMatrixLumped(self, matrix):
    """
    Calculate lumped matrix associated with integration.
    """

    (nrows, ncols) = matrix.shape
    lumped = numpy.zeros( (ncols), dtype=numpy.float64)
    
    for iVertex in xrange(self.numCorners):
      i = numpy.asarray(range(self.numCorners))
      for iDim in xrange(self.fiberDim):
        iR = iVertex * self.fiberDim + iDim
        indices = self.fiberDim * i + iDim
        lumped[iR] = sum(matrix[iR, indices])
    return lumped


# End of file 
