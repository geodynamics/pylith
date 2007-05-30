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

## @file unittests/libtests/feassemble/data/Quadrature3DLinear.odb
##
## @brief Python container holding quadrature information for a 3-D
## linear finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------
def N0(p):
  return 1.0-p[0]-p[1]-p[2]

def N0p(p):
  return -1.0

def N0q(p):
  return -1.0

def N0r(p):
  return -1.0

def N1(p):
  return p[0]

def N1p(p):
  return 1.0

def N1q(p):
  return 0.0

def N1r(p):
  return 0.0

def N2(p):
  return p[1]

def N2p(p):
  return 0.0

def N2q(p):
  return 1.0

def N2r(p):
  return 0.0

def N3(p):
  return p[2]

def N3p(p):
  return 0.0

def N3q(p):
  return 0.0

def N3r(p):
  return 1.0


# ----------------------------------------------------------------------

# Quadrature3DLinear class
class Quadrature3DLinear(Component):
  """
  Python container holding quadrature information for a 3-D linear
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="quadrature3dlinear"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="quadrature")
    
    self.quadPtsRef = numpy.array( [[1.0/4.0, 1.0/4.0, 1.0/4.0]],
                                   dtype=numpy.float64)
    self.quadWts = numpy.array([1.0/6.0], dtype=numpy.float64)
    self.numBasis = 4
    self.numQuadPts = 1
    self.spaceDim = 3
    self.cellDim = 3
    return
  

  def calculateBasis(self):
    """
    Calculate basis functions and their derivatives at quadrature points.
    """

    basis = numpy.zeros( (self.numQuadPts, self.numBasis),
                         dtype=numpy.float64)
    basisDeriv = numpy.zeros( (self.numQuadPts, self.numBasis, self.cellDim),
                              dtype=numpy.float64)

    iQuad = 0
    for q in self.quadPtsRef:
      # Basis functions at quadrature points
      basisQ = numpy.array([N0(q), N1(q), N2(q), N3(q)], dtype=numpy.float64)
      basis[iQuad] = basisQ.reshape( (self.numBasis,) )

      # Derivatives of basis functions at quadrature points
      derivQ = numpy.array([[N0p(q), N0q(q), N0r(q)],
                            [N1p(q), N1q(q), N1r(q)],
                            [N2p(q), N2q(q), N2r(q)],
                            [N3p(q), N3q(q), N3r(q)]],
                           dtype=numpy.float64)      
      basisDeriv[iQuad] = derivQ.reshape((self.numBasis, self.cellDim))

      iQuad += 1
    return (basis, basisDeriv)
    

# FACTORIES ////////////////////////////////////////////////////////////
def quadrature():
  """
  Factory for Quadrature3DLinear.
  """
  return Quadrature3DLinear()


# End of file 
