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

## @file unittests/libtests/feassemble/data/Quadrature3DQuadratic.py

## @brief Python application for generating C++ data files for testing
## C++ Quadrature3D object w/quadratic basis functions.

from QuadratureApp import QuadratureApp

import numpy

# ----------------------------------------------------------------------
def N0(p):
  return (1.0-p[0]-p[1]-p[2])*(1.0-2.0*p[0]-2.0*p[1]-2.0*p[2])

def N0p(p):
  return 4.0*p[0]+4.0*p[1]+4.0*p[2]-3.0

def N0q(p):
  return 4.0*p[0]+4.0*p[1]+4.0*p[2]-3.0

def N0r(p):
  return 4.0*p[0]+4.0*p[1]+4.0*p[2]-3.0

def N1(p):
  return p[0]*(2.0*p[0]-1.0)

def N1p(p):
  return 4.0*p[0]-1

def N1q(p):
  return 0.0

def N1r(p):
  return 0.0

def N2(p):
  return p[1]*(2.0*p[1]-1.0)

def N2p(p):
  return 0.0

def N2q(p):
  return 4.0*p[1]-1.0

def N2r(p):
  return 0.0

def N3(p):
  return p[2]*(2.0*p[2]-1.0)

def N3p(p):
  return 0.0

def N3q(p):
  return 0.0

def N3r(p):
  return 4.0*p[2]-1.0

def N4(p):
  return 4.0*p[0]*(1.0-p[0]-p[1]-p[2])

def N4p(p):
  return -8.0*p[0]+4.0*(1.0-p[1]-p[2])

def N4q(p):
  return -4.0*p[0]

def N4r(p):
  return -4.0*p[0]

def N5(p):
  return 4.0*p[1]*(1.0-p[0]-p[1]-p[2])

def N5p(p):
  return -4.0*p[1]

def N5q(p):
  return -8.0*p[1]+4.0*(1.0-p[0]-p[2])

def N5r(p):
  return -4.0*p[1]

def N6(p):
  return 4.0*p[2]*(1.0-p[0]-p[1]-p[2])

def N6p(p):
  return -4.0*p[2]

def N6q(p):
  return -4.0*p[2]

def N6r(p):
  return -8.0*p[2]+4.0*(1.0-p[0]-p[1])

def N7(p):
  return 4.0*p[0]*p[1]

def N7p(p):
  return 4.0*p[1]

def N7q(p):
  return 4.0*p[0]

def N7r(p):
  return 0.0

def N8(p):
  return 4.0*p[1]*p[2]

def N8p(p):
  return 0.0

def N8q(p):
  return 4.0*p[2]

def N8r(p):
  return 4.0*p[1]

def N9(p):
  return 4.0*p[0]*p[2]

def N9p(p):
  return 4.0*p[2]

def N9q(p):
  return 0.0

def N9r(p):
  return 4.0*p[0]

# ----------------------------------------------------------------------

# Quadrature3DQuadratic class
class Quadrature3DQuadratic(QuadratureApp):
  """
  Python application for generating C++ data files for testing C++
  Quadrature3D object w/quadratic basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature3dquadratic"):
    """
    Constructor.
    """
    QuadratureApp.__init__(self, name)
    
    self.numVertices = 10
    self.spaceDim = 3
    self.numCells = 1
    self.cellDim = 3
    self.numCorners = 10
    self.numQuadPts = 4

    # These are just approximate points used to test the quadrature routine
    self.quadPtsRef = numpy.array( [[1.0/12.0, 1.0/12.0, 1.0/12.0],
                                    [3.0/4.0, 1.0/12.0, 1.0/12.0],
                                    [1.0/12.0, 3.0/4.0, 1.0/12.0],
                                    [1.0/12.0, 1.0/12.0, 3.0/4.0]],
                                   dtype=numpy.Float64)
    self.quadWts = numpy.array([1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0],
                               dtype=numpy.Float64)
    self.vertices = numpy.array( [[0.0, 0.0, 0.0],
                                  [1.0, 0.0, 0.0],
                                  [0.0, 1.0, 0.0],
                                  [0.0, 0.0, 1.0],
                                  [0.5, 0.0, 0.0],
                                  [0.0, 0.5, 0.0],
                                  [0.0, 0.0, 0.5],
                                  [0.5, 0.5, 0.0],
                                  [0.0, 0.5, 0.5],
                                  [0.5, 0.0, 0.5]],
                                 dtype=numpy.Float64)
    
    self.vertices = numpy.array( [[-0.5, -2.0, -1.0],
                                  [ 2.0, -2.0, -0.5],
                                  [ 1.0,  1.0,  0.0],
                                  [ 0.2,  0.5,  2.0],
                                  [ 0.7, -2.1, -0.8],
                                  [ 0.3, -0.5, -0.5],
                                  [-0.2, -0.8,  0.5],
                                  [ 1.5, -0.6, -0.2],
                                  [ 0.6,  0.8,  0.9],
                                  [ 1.1, -0.8,  0.7]],
                                 dtype=numpy.Float64)
    self.cells = numpy.array( [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]],
                              dtype=numpy.Int32)
    
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
      basis = numpy.array([N0(q), N1(q), N2(q), N3(q), N4(q),
                           N5(q), N6(q), N7(q), N8(q), N9(q)],
                          dtype=numpy.Float64)
      self.basis[iQuad] = basis.reshape( (self.numCorners,) )

      # Derivatives of basis functions at quadrature points
      deriv = numpy.array([[N0p(q), N0q(q), N0r(q)],
                           [N1p(q), N1q(q), N1r(q)],
                           [N2p(q), N2q(q), N2r(q)],
                           [N3p(q), N3q(q), N3r(q)],
                           [N4p(q), N4q(q), N4r(q)],
                           [N5p(q), N5q(q), N5r(q)],
                           [N6p(q), N6q(q), N6r(q)],
                           [N7p(q), N7q(q), N7r(q)],
                           [N8p(q), N8q(q), N8r(q)],
                           [N9p(q), N9q(q), N9r(q)]],
                          dtype=numpy.Float64)      
      self.basisDeriv[iQuad] = deriv.reshape((self.numCorners, self.cellDim))

      iQuad += 1
    return
    

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = Quadrature3DQuadratic()
  app.run()


# End of file 
