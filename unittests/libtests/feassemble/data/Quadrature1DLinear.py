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

## @file unittests/libtests/feassemble/data/Quadrature1DLinear.py

## @brief Python application for generating C++ data files for testing
## C++ Quadrature1D object w/linear basis functions.

from QuadratureApp import QuadratureApp

import numpy

# ----------------------------------------------------------------------
def N0(p):
  return 0.5*(1.0-p)

def N0p(p):
  return -0.5

def N1(p):
  return 0.5*(1.0+p)

def N1p(p):
  return +0.5

# ----------------------------------------------------------------------

# Quadrature1DLinear class
class Quadrature1DLinear(QuadratureApp):
  """
  Python application for generating C++ data files for testing C++
  Quadrature1D object w/linear basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="quadrature1dlinear"):
    """
    Constructor.
    """
    QuadratureApp.__init__(self, name)
    
    self.numVertices = 2
    self.spaceDim = 1
    self.numCells = 1
    self.cellDim = 1
    self.numCorners = 2
    self.numQuadPts = 1
    
    self.quadPtsRef = numpy.array( [[0.0]], dtype=numpy.Float64)
    self.quadWts = numpy.array([2.0], dtype=numpy.Float64)
    self.vertices = numpy.array( [[-0.25], [2.0]], dtype=numpy.Float64)
    self.cells = numpy.array( [[0, 1]], dtype=numpy.Int32)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////
  
  def _calculateBasis(self):
    """
    Calculate basis functions and their derivatives at quadrature points.
    """

    self.basis = numpy.zeros( (self.numQuadPts, self.numCorners),
                              dtype=numpy.Float64)
    self.basisDeriv = numpy.zeros( (self.numQuadPts,
                                    self.numCorners, self.cellDim),
                                   dtype=numpy.Float64)

    iQuad = 0
    for q in self.quadPtsRef:
      # Basis functions at quadrature points
      basis = numpy.array([N0(q), N1(q)], dtype=numpy.Float64)
      self.basis[iQuad,:] = basis.reshape( (self.numCorners,) )

      # Derivatives of basis functions at quadrature points
      deriv = numpy.array([[N0p(q)], [N1p(q)]], dtype=numpy.Float64)      
      self.basisDeriv[iQuad,:] = deriv.reshape((self.numCorners, self.cellDim))

      iQuad += 1
    return
    

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = Quadrature1DLinear()
  app.run()


# End of file 
