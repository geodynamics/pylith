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

## @file unittests/libtests/feassemble/data/Quadrature1DLinear.odb
##
## @brief Python container holding quadrature information for a 1-D
## linear finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

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
class Quadrature1DLinear(Component):
  """
  Python container holding quadrature information for a 1-D linear
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="quadrature1dlinear"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="quadrature")
    
    self.quadPtsRef = numpy.array( [[0.0]], dtype=numpy.float64)
    self.quadWts = numpy.array([2.0], dtype=numpy.float64)
    self.numBasis = 2
    self.numQuadPts = 1
    self.spaceDim = 1
    self.cellDim = 1
    return
  

  def calculateBasis(self):
    """
    Calculate basis functions and their derivatives at quadrature points.
    """

    basis = numpy.zeros( (self.numQuadPts, self.numBasis), dtype=numpy.float64)
    basisDeriv = numpy.zeros( (self.numQuadPts, self.numBasis, self.cellDim),
                              dtype=numpy.float64)

    iQuad = 0
    for q in self.quadPtsRef:
      # Basis functions at quadrature points
      basisQ = numpy.array([N0(q), N1(q)], dtype=numpy.float64)
      basis[iQuad,:] = basisQ.reshape( (self.numBasis,) )

      # Derivatives of basis functions at quadrature points
      derivQ = numpy.array([[N0p(q)], [N1p(q)]], dtype=numpy.float64)      
      basisDeriv[iQuad,:] = derivQ.reshape((self.numBasis, self.cellDim))

      iQuad += 1
    return (basis, basisDeriv)


# FACTORIES ////////////////////////////////////////////////////////////
def quadrature():
  """
  Factory for Quadrature1DLinear.
  """
  return Quadrature1DLinear()


# End of file 
