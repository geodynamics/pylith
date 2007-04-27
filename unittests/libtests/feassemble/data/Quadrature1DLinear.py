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

def verticesRef():
  return [-1.0, 1.0]

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
    self.numBasis = 2
    self.numQuadPts = 1
    
    self.quadPtsRef = numpy.array( [[0.0]], dtype=numpy.float64)
    self.quadWts = numpy.array([2.0], dtype=numpy.float64)
    self.vertices = numpy.array( [[-0.25], [2.0]], dtype=numpy.float64)
    self.cells = numpy.array( [[0, 1]], dtype=numpy.int32)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////
  
  def calculateBasisVert(self):
    """
    Calculate basis functions and their derivatives at reference vertices.
    """

    self.basisVert = numpy.zeros( (self.numBasis, self.numBasis),
                                  dtype=numpy.float64)
    self.basisDerivVert = numpy.zeros( (self.numBasis,
                                        self.numBasis, self.cellDim),
                                       dtype=numpy.float64)

    iVertex = 0
    for v in verticesRef():
      # Basis functions at vertices
      basis = numpy.array([N0(v), N1(v)], dtype=numpy.float64)
      self.basisVert[iVertex,:] = basis.reshape( (self.numBasis,) )

      # Derivatives of basis functions at quadrature points
      deriv = numpy.array([[N0p(v)], [N1p(v)]], dtype=numpy.float64)      
      self.basisDerivVert[iVertex,:] = deriv.reshape((self.numBasis,
                                                      self.cellDim))

      iVertex += 1
    return
    

  def calculateBasisQuad(self):
    """
    Calculate basis functions and their derivatives at quadrature points.
    """

    self.basisQuad = numpy.zeros( (self.numQuadPts, self.numBasis),
                                  dtype=numpy.float64)
    self.basisDerivQuad = numpy.zeros( (self.numQuadPts,
                                        self.numBasis, self.cellDim),
                                       dtype=numpy.float64)

    iQuad = 0
    for q in self.quadPtsRef:
      # Basis functions at quadrature points
      basis = numpy.array([N0(q), N1(q)], dtype=numpy.float64)
      self.basisQuad[iQuad,:] = basis.reshape( (self.numBasis,) )

      # Derivatives of basis functions at quadrature points
      deriv = numpy.array([[N0p(q)], [N1p(q)]], dtype=numpy.float64)      
      self.basisDerivQuad[iQuad,:] = deriv.reshape((self.numBasis,
                                                    self.cellDim))

      iQuad += 1
    return
    

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = Quadrature1DLinear()
  app.run()


# End of file 
