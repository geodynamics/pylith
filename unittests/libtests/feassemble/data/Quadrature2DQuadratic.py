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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/Quadrature2DQuadratic.odb
##
## @brief Python container holding quadrature information for a 2-D
## quadratic finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------
def N0(p):
  return 0.5*(1.0+p[0]+p[1])*(p[0]+p[1])

def N0p(p):
  return 0.5+p[0]+p[1]

def N0q(p):
  return 0.5+p[0]+p[1]

def N1(p):
  return 0.5*p[0]*(p[0]+1.0)

def N1p(p):
  return 0.5+p[0]

def N1q(p):
  return 0.0

def N2(p):
  return 0.5*p[1]*(p[1]+1.0)

def N2p(p):
  return 0.0

def N2q(p):
  return 0.5+p[1]

def N3(p):
  return (1.0+p[0])*(1.0+p[1])

def N3p(p):
  return 1.0+p[1]

def N3q(p):
  return 1.0+p[0]

def N4(p):
  return -(p[0]+p[1])*(1.0+p[1])

def N4p(p):
  return -(1.0+p[1])

def N4q(p):
  return -1.0-p[0]-2.0*p[1]

def N5(p):
  return -(p[0]+p[1])*(1.0+p[0])

def N5p(p):
  return -1.0-p[1]-2.0*p[0]

def N5q(p):
  return -(1.0+p[0])

# ----------------------------------------------------------------------

# Quadrature2DQuadratic class
class Quadrature2DQuadratic(Component):
  """
  Python container holding quadrature information for a 1-D quadratic
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature2dquadratic"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="quadrature")
    
    # These are just approximate points used to test the quadrature routine
    self.quadPtsRef = numpy.array( [[-0.75,-0.75],
                                    [0.75,-0.75],
                                    [-0.75,0.75],
                                    [0,-0.75],
                                    [-0.75,0],
                                    [0.25,0.25]],
                                   dtype=numpy.float64)
    self.quadWts = numpy.array([1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0],
                               dtype=numpy.float64)
    #self.quadPtsRef = numpy.array( [[-0.64288254347276719, -0.68989794855663567],
    #                                [-0.84993777955478378, 0.28989794855663559],
    #                                [0.33278049202940285, -0.68989794855663567],
    #                                [-0.43996016900185175, 0.28989794855663559]],
    #                               dtype=numpy.float64)
    #self.quadWts = numpy.array([0.63608276,  0.36391724,  0.63608276,  0.36391724],
    #                           dtype=numpy.float64)
    self.numBasis = 6
    self.numQuadPts = 6
    self.spaceDim = 2
    self.cellDim = 2
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
      basisQ = numpy.array([N0(q), N1(q), N2(q), N3(q), N4(q), N5(q)],
                          dtype=numpy.float64)
      basis[iQuad] = basisQ.reshape( (self.numBasis,) )
      
      # Derivatives of basis functions at quadrature points
      derivQ = numpy.array([[N0p(q), N0q(q)],
                           [N1p(q), N1q(q)],
                           [N2p(q), N2q(q)],
                           [N3p(q), N3q(q)],
                           [N4p(q), N4q(q)],
                           [N5p(q), N5q(q)]],
                          dtype=numpy.float64)      
      basisDeriv[iQuad] = derivQ.reshape((self.numBasis, self.cellDim))
      
      iQuad += 1
    return (basis, basisDeriv)
    

# FACTORIES ////////////////////////////////////////////////////////////
def quadrature():
  """
  Factory for Quadrature2DQuadratic.
  """
  return Quadrature2DQuadratic()


# End of file 
