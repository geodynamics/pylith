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

## @file unittests/libtests/feassemble/data/Quadrature1DQuadratic.py

## @brief Python application for generating C++ data files for testing
## C++ Quadrature1D object w/quadratic basis functions.

from QuadratureApp import QuadratureApp

import numpy

# ----------------------------------------------------------------------
def N0(p):
  return -0.5*p*(1.0-p)

def N0p(p):
  return -0.5*(1.0-p) + 0.5*p

def N1(p):
  return (1.0-p**2)

def N1p(p):
  return -2.0*p

def N2(p):
  return 0.5*p*(1.0+p)

def N2p(p):
  return +0.5*(1.0+p) + 0.5*p

# ----------------------------------------------------------------------

# Quadrature1DQuadratic class
class Quadrature1DQuadratic(QuadratureApp):
  """
  Python application for generating C++ data files for testing C++
  Quadrature1D object w/quadratic basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature1dquadratic"):
    """
    Constructor.
    """
    QuadratureApp.__init__(self, name)
    
    self.numVertices = 3
    self.spaceDim = 1
    self.numCells = 1
    self.cellDim = 1
    self.numCorners = 3
    self.numQuadPts = 2
    
    self.quadPtsRef = numpy.array( [[-1.0/3**0.5],
                                    [+1.0/3**0.5]],
                                   dtype=numpy.Float64)
    self.quadWts = numpy.array([1.0, 1.0], dtype=numpy.Float64)
    self.vertices = numpy.array( [[-0.25], [0.875], [2.0]],
                                 dtype=numpy.Float64)
    self.cells = numpy.array( [[0, 1, 2]], dtype=numpy.Int32)
    
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////
  
  def _calculate(self):
    """
    Calculate basis functions, derivatives, and Jacobian information
    at quadrature points.
    """

    # Basis functions at quadrature points
    n0 = N0(self.quadPtsRef)
    n1 = N1(self.quadPtsRef)
    n2 = N2(self.quadPtsRef)
    basis = numpy.array([n0, n1, n2]).transpose()
    self.basis = basis.reshape( (self.numQuadPts, self.numVertices) )

    # Derivatives of basis functions at quadrature points
    n0p = N0p(self.quadPtsRef)
    n1p = N1p(self.quadPtsRef)
    n2p = N2p(self.quadPtsRef)
    deriv = numpy.array([n0p, n1p, n2p]).transpose()
    self.basisDeriv = deriv.reshape( (self.numQuadPts, self.numVertices) )

    # Jacobian at quadrature points
    self.jacobian = numpy.dot(self.basisDeriv, self.vertices)

    # Determinant of Jacobian at quadrature points
    self.jacobianDet = self.jacobian

    # Inverse of Jacobian at quadrature points
    self.jacobianInv = 1.0 / self.jacobian

    # Quadrature points in cell
    self.quadPts = numpy.dot(self.basis, self.vertices)
    return
      

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = Quadrature1DQuadratic()
  app.run()


# End of file 
