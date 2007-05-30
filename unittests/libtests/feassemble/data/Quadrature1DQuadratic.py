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

## @file unittests/libtests/feassemble/data/Quadrature1DQuadratic.odb
##
## @brief Python application for generating C++ data files for testing
## C++ Quadrature1D object w/quadratic basis functions.

from pyre.components.Component import Component

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
class Quadrature1DQuadratic(Component):
  """
  Python application for generating C++ data files for testing C++
  Quadrature1D object w/quadratic basis functions.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature1dquadratic"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="quadrature")
    
    self.quadPtsRef = numpy.array( [[-1.0/3**0.5],
                                    [+1.0/3**0.5]],
                                   dtype=numpy.float64)
    self.quadWts = numpy.array([1.0, 1.0], dtype=numpy.float64)
    self.numBasis = 3
    self.numQuadPts = 2
    self.spaceDim = 1
    self.cellDim = 1    
    return


  def calculateBasis(self):
    """
    Calculate basis functions and their derivatives at quadrature points.
    """

    basis = numpy.zeros( (self.numQuadPts, self.numBasis),
                         dtype=numpy.float64)
    basisDeriv = numpy.zeros( (self.numQuadPts,
                               self.numBasis, self.cellDim),
                              dtype=numpy.float64)

    iQuad = 0
    for q in self.quadPtsRef:
      # Basis functions at quadrature points
      basisQ = numpy.array([N0(q), N1(q), N2(q)], dtype=numpy.float64)
      basis[iQuad] = basisQ.reshape( (self.numBasis,) )

      # Derivatives of basis functions at quadrature points
      derivQ = numpy.array([[N0p(q)], [N1p(q)], [N2p(q)]],
                          dtype=numpy.float64)      
      basisDeriv[iQuad] = derivQ.reshape((self.numBasis, self.cellDim))

      iQuad += 1
    return (basis, basisDeriv)
    

# FACTORIES ////////////////////////////////////////////////////////////
def quadrature():
  """
  Factory for Quadrature1DLinear.
  """
  return Quadrature1DQuadratic()


# End of file 
