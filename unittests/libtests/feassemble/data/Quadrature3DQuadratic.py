#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/Quadrature3DQuadratic.odb
##
## @brief Python container holding quadrature information for a 3-D
## quadratic finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------
def N0(p):
  return 0.5 * (1.0+p[0]+p[1]+p[2])*(2.0+p[0]+p[1]+p[2])

def N0p(p):
  return 1.5+p[0]+p[1]+p[2]

def N0q(p):
  return 1.5+p[0]+p[1]+p[2]

def N0r(p):
  return 1.5+p[0]+p[1]+p[2]

def N1(p):
  return 0.5*(1.0+p[0])*p[0]

def N1p(p):
  return 0.5+p[0]

def N1q(p):
  return 0.0

def N1r(p):
  return 0.0

def N2(p):
  return 0.5*(1.0+p[1])*p[1]

def N2p(p):
  return 0.0

def N2q(p):
  return 0.5+p[1]

def N2r(p):
  return 0.0

def N3(p):
  return 0.5*(1.0+p[2])*p[2]

def N3p(p):
  return 0.0

def N3q(p):
  return 0.0

def N3r(p):
  return 0.5+p[2]

def N4(p):
  return (1.0+p[0])*(1.0+p[1])

def N4p(p):
  return 1.0+p[1]

def N4q(p):
  return 1.0+p[0]

def N4r(p):
  return 0.0

def N5(p):
  return -(1.0+p[1])*(1.0+p[0]+p[1]+p[2])

def N5p(p):
  return -(1.0+p[1])

def N5q(p):
  return -2.0-p[0]-2.0*p[1]-p[2]

def N5r(p):
  return -(1.0+p[1])

def N6(p):
  return -(1.0+p[0])*(1.0+p[0]+p[1]+p[2])

def N6p(p):
  return -2.0-2.0*p[0]-p[1]-p[2]

def N6q(p):
  return -(1.0+p[0])

def N6r(p):
  return -(1.0+p[0])

def N7(p):
  return -(1.0+p[2])*(1.0+p[0]+p[1]+p[2])

def N7p(p):
  return -(1.0+p[2])

def N7q(p):
  return -(1.0+p[2])

def N7r(p):
  return -2.0-p[0]-p[1]-2.0*p[2]

def N8(p):
  return (1.0+p[0])*(1.0+p[2])

def N8p(p):
  return 1.0 + p[2]

def N8q(p):
  return 0.0

def N8r(p):
  return 1.0 + p[0]

def N9(p):
  return (1.0+p[1])*(1.0+p[2])

def N9p(p):
  return 0.0

def N9q(p):
  return 1.0 + p[2]

def N9r(p):
  return 1.0 + p[1]


# ----------------------------------------------------------------------

# Quadrature3DQuadratic class
class Quadrature3DQuadratic(Component):
  """
  Python container holding quadrature information for a 3-D quadratic
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature3dquadratic"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="quadrature")
    

    # These are just approximate points used to test the quadrature routine
    self.quadPtsRef = numpy.array( [[-1.0+1.0/5.0, -1.0+1.0/5.0, -1.0+1.0/5.0],
                                    [-1.0+3.0/2.0, -1.0+1.0/5.0, -1.0+1.0/5.0],
                                    [-1.0+1.0/5.0, -1.0+3.0/2.0, -1.0+1.0/5.0],
                                    [-1.0+1.0/5.0, -1.0+1.0/5.0, -1.0+3.0/2.0]],
                                   dtype=numpy.float64)
    self.quadWts = numpy.array([1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0],
                               dtype=numpy.float64)
    self.numBasis = 10
    self.numQuadPts = 4
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
      basisQ = numpy.array([N0(q), N1(q), N2(q), N3(q), N4(q),
                            N5(q), N6(q), N7(q), N8(q), N9(q)],
                           dtype=numpy.float64)
      basis[iQuad] = basisQ.reshape( (self.numBasis,) )

      # Derivatives of basis functions at quadrature points
      derivQ = numpy.array([[N0p(q), N0q(q), N0r(q)],
                            [N1p(q), N1q(q), N1r(q)],
                            [N2p(q), N2q(q), N2r(q)],
                            [N3p(q), N3q(q), N3r(q)],
                            [N4p(q), N4q(q), N4r(q)],
                            [N5p(q), N5q(q), N5r(q)],
                            [N6p(q), N6q(q), N6r(q)],
                            [N7p(q), N7q(q), N7r(q)],
                            [N8p(q), N8q(q), N8r(q)],
                            [N9p(q), N9q(q), N9r(q)]],
                           dtype=numpy.float64)      
      basisDeriv[iQuad] = derivQ.reshape((self.numBasis, self.cellDim))

      iQuad += 1
    return (basis, basisDeriv)
    

# FACTORIES ////////////////////////////////////////////////////////////
def quadrature():
  """
  Factory for Quadrature3DQuadratic.
  """
  return Quadrature3DQuadratic()


# End of file 
