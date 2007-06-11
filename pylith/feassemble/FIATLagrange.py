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

## @file pylith/feassemble/FIATLagrange.py
##
## @brief Python object for managing basis functions and quadrature
## rules of a Lagrange reference finite-element cell using FIAT.
##
## The basis functions are constructed from the tensor product of 1-D
## Lagrance reference cells.
##
## Factory: reference_cell.

from ReferenceCell import ReferenceCell

import numpy

def validateDimension(dim):
  if dim < 1 or dim > 3:
    raise ValueError("Dimension of Lagrange element must be 1, 2, or 3.")
  return dim

# FIATLagrange class
class FIATLagrange(ReferenceCell):
  """
  Python object for managing basis functions and quadrature rules of a
  Lagrange reference finite-element cell using FIAT.

  Factory: reference_cell.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ReferenceCell.Inventory):
    """Python object for managing FIATLagrange facilities and properties."""

    ## @class Inventory
    ## Python object for managing FIATLagrange facilities and properties.
    ##
    ## \b Properties
    ## @li \b dimension Dimension of finite-element cell
    ## @li \b degree Degree of finite-element cell 
    ## @li \b quad_order Order of quadrature rule
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    dimension = pyre.inventory.int("dimension", default=3,
                                   validator=validateDimension)
    dimension.meta['tip'] = "Dimension of finite-element cell."

    degree = pyre.inventory.int("degree", default=1)
    degree.meta['tip'] = "Degree of finite-element cell."

    order = pyre.inventory.int("quad_order", default=-1)
    order.meta['tip'] = "Order of quadrature rule."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fiatlagrange"):
    """
    Constructor.
    """
    ReferenceCell.__init__(self, name)
    return


  def initialize(self, spaceDim):
    """
    Initialize reference finite-element cell from a tensor product of
    1-D Lagrange elements.
    """
    self._setupGeometry(spaceDim)

    if  self.cellDim > 0:
      quadrature = self._setupQuadrature()
      element    = self._setupElement()
      dim        = self.cellDim
    
      # Get coordinates of vertices (dual basis)
      vertices = numpy.array(self._setupVertices(element))

      # Evaluate basis functions at quadrature points
      quadpts = numpy.array(quadrature.get_points())
      quadwts = numpy.array(quadrature.get_weights())
      numQuadPts = len(quadpts)
      basis = numpy.array(element.function_space().tabulate(quadrature.get_points())).transpose()
      numBasisFns = len(element.function_space())

      # Evaluate derivatives of basis functions at quadrature points
      basisDeriv = numpy.array([element.function_space().deriv_all(d).tabulate(quadrature.get_points()) \
                                for d in range(1)]).transpose()

      self.numQuadPts = numQuadPts**dim
      self.numCorners = numBasisFns**dim

      if dim == 1:
        self.vertices = numpy.array(vertices)
        self.quadPts = quadpts
        self.quadWts = quadwts
        self.basis = basis
        self.basisDeriv = basisDeriv
      else:
        if dim == 2:
          self.vertices = numpy.zeros((self.numCorners, dim))
          n = 0
          # Bottom
          for r in range(0, numBasisFns-1):
            self.vertices[n][0] = vertices[r]
            self.vertices[n][1] = vertices[0]
            n += 1
          # Right
          for q in range(0, numBasisFns-1):
            self.vertices[n][0] = vertices[numBasisFns-1]
            self.vertices[n][1] = vertices[q]
            n += 1
          # Top
          for r in range(numBasisFns-1, 0, -1):
            self.vertices[n][0] = vertices[r]
            self.vertices[n][1] = vertices[numBasisFns-1]
            n += 1
          # Left
          for q in range(numBasisFns-1, 0, -1):
            self.vertices[n][0] = vertices[0]
            self.vertices[n][1] = vertices[q]
            n += 1
          # Interior
          for q in range(1, numBasisFns-1):
            for r in range(1, numBasisFns-1):
              self.vertices[n][0] = vertices[r]
              self.vertices[n][1] = vertices[q]
              n += 1
          if not n == self.numCorners: raise RuntimeError('Invalid 2D function tabulation')
        
          self.quadPts = numpy.zeros((numQuadPts*numQuadPts, dim))
          self.quadWts = numpy.zeros((numQuadPts*numQuadPts,))
          self.basis = numpy.zeros((numQuadPts*numQuadPts,
                                         numBasisFns*numBasisFns))
          self.basisDeriv = numpy.zeros((numQuadPts*numQuadPts,
                                         numBasisFns*numBasisFns, dim))
          n = 0
          # Bottom
          for r in range(0, numQuadPts-1):
            self.quadPts[n][0] = quadpts[r]
            self.quadPts[n][1] = quadpts[0]
            self.quadWts[n]    = quadwts[r]*quadwts[0]
            m = 0
            # Bottom
            for g in range(0, numBasisFns-1):
              self.basis[n][m] = basis[r][g]*basis[0][0]
              self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[0][0]
              self.basisDeriv[n][m][1] = basis[r][g]*basisDeriv[0][0][0]
              m += 1
            # Right
            for f in range(0, numBasisFns-1):
              self.basis[n][m] = basis[r][numBasisFns-1]*basis[0][f]
              self.basisDeriv[n][m][0] = basisDeriv[r][numBasisFns-1][0]*basis[0][f]
              self.basisDeriv[n][m][1] = basis[r][numBasisFns-1]*basisDeriv[0][f][0]
              m += 1
            # Top
            for g in range(numBasisFns-1, 0, -1):
              self.basis[n][m] = basis[r][g]*basis[0][numBasisFns-1]
              self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[0][numBasisFns-1]
              self.basisDeriv[n][m][1] = basis[r][g]*basisDeriv[0][numBasisFns-1][0]
              m += 1
            # Left
            for f in range(numBasisFns-1, 0, -1):
              self.basis[n][m] = basis[r][0]*basis[0][f]
              self.basisDeriv[n][m][0] = basisDeriv[r][0][0]*basis[0][f]
              self.basisDeriv[n][m][1] = basis[r][0]*basisDeriv[0][f][0]
              m += 1
            # Interior
            for f in range(1, numBasisFns-1):
              for g in range(1, numBasisFns-1):
                self.basis[n][m] = basis[r][g]*basis[0][f]
                self.basisDeriv[0][r][f][g][0] = basisDeriv[r][g][0]*basis[0][f]
                self.basisDeriv[0][r][f][g][1] = basis[r][g]*basisDeriv[0][f][0]
                m += 1
            if not m == self.numCorners: raise RuntimeError('Invalid 2D function tabulation')
            n += 1
          # Right
          for q in range(0, numQuadPts-1):
            self.quadPts[n][0] = quadpts[numQuadPts-1]
            self.quadPts[n][1] = quadpts[q]
            self.quadWts[n]    = quadwts[numQuadPts-1]*quadwts[q]
            m = 0
            # Bottom
            for g in range(0, numBasisFns-1):
              self.basis[n][m] = basis[numQuadPts-1][g]*basis[q][0]
              self.basisDeriv[n][m][0] = basisDeriv[numQuadPts-1][g][0]*basis[q][0]
              self.basisDeriv[n][m][1] = basis[numQuadPts-1][g]*basisDeriv[q][0][0]
              m += 1
            # Right
            for f in range(0, numBasisFns-1):
              self.basis[n][m] = basis[numQuadPts-1][numBasisFns-1]*basis[q][f]
              self.basisDeriv[n][m][0] = basisDeriv[numQuadPts-1][numBasisFns-1][0]*basis[q][f]
              self.basisDeriv[n][m][1] = basis[numQuadPts-1][numBasisFns-1]*basisDeriv[q][f][0]
              m += 1
            # Top
            for g in range(numBasisFns-1, 0, -1):
              self.basis[n][m] = basis[numQuadPts-1][g]*basis[q][numBasisFns-1]
              self.basisDeriv[n][m][0] = basisDeriv[numQuadPts-1][g][0]*basis[q][numBasisFns-1]
              self.basisDeriv[n][m][1] = basis[numQuadPts-1][g]*basisDeriv[q][numBasisFns-1][0]
              m += 1
            # Left
            for f in range(numBasisFns-1, 0, -1):
              self.basis[n][m] = basis[numQuadPts-1][0]*basis[q][f]
              self.basisDeriv[n][m][0] = basisDeriv[numQuadPts-1][0][0]*basis[q][f]
              self.basisDeriv[n][m][1] = basis[numQuadPts-1][0]*basisDeriv[q][f][0]
              m += 1
            # Interior
            for f in range(1, numBasisFns-1):
              for g in range(1, numBasisFns-1):
                self.basis[n][m] = basis[numQuadPts-1][g]*basis[0][f]
                self.basisDeriv[q][numQuadPts-1][f][g][0] = basisDeriv[numQuadPts-1][g][0]*basis[q][f]
                self.basisDeriv[q][numQuadPts-1][f][g][1] = basis[numQuadPts-1][g]*basisDeriv[q][f][0]
                m += 1
            if not m == self.numCorners: raise RuntimeError('Invalid 2D function tabulation')
            n += 1
          # Top
          for r in range(numQuadPts-1, 0, -1):
            self.quadPts[n][0] = quadpts[r]
            self.quadPts[n][1] = quadpts[numQuadPts-1]
            self.quadWts[n]    = quadwts[r]*quadwts[numQuadPts-1]
            m = 0
            # Bottom
            for g in range(0, numBasisFns-1):
              self.basis[n][m] = basis[r][g]*basis[numQuadPts-1][0]
              self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[numQuadPts-1][0]
              self.basisDeriv[n][m][1] = basis[r][g]*basisDeriv[numQuadPts-1][0][0]
              m += 1
            # Right
            for f in range(0, numBasisFns-1):
              self.basis[n][m] = basis[r][numBasisFns-1]*basis[numQuadPts-1][f]
              self.basisDeriv[n][m][0] = basisDeriv[r][numBasisFns-1][0]*basis[numQuadPts-1][f]
              self.basisDeriv[n][m][1] = basis[r][numBasisFns-1]*basisDeriv[numQuadPts-1][f][0]
              m += 1
            # Top
            for g in range(numBasisFns-1, 0, -1):
              self.basis[n][m] = basis[r][g]*basis[numQuadPts-1][numBasisFns-1]
              self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[numQuadPts-1][numBasisFns-1]
              self.basisDeriv[n][m][1] = basis[r][g]*basisDeriv[numQuadPts-1][numBasisFns-1][0]
              m += 1
            # Left
            for f in range(numBasisFns-1, 0, -1):
              self.basis[n][m] = basis[r][0]*basis[numQuadPts-1][f]
              self.basisDeriv[n][m][0] = basisDeriv[r][0][0]*basis[numQuadPts-1][f]
              self.basisDeriv[n][m][1] = basis[r][0]*basisDeriv[numQuadPts-1][f][0]
              m += 1
            # Interior
            for f in range(1, numBasisFns-1):
              for g in range(1, numBasisFns-1):
                self.basis[n][m] = basis[r][g]*basis[numQuadPts-1][f]
                self.basisDeriv[numQuadPts-1][r][f][g][0] = basisDeriv[r][g][0]*basis[numQuadPts-1][f]
                self.basisDeriv[numQuadPts-1][r][f][g][1] = basis[r][g]*basisDeriv[numQuadPts-1][f][0]
                m += 1
            if not m == self.numCorners: raise RuntimeError('Invalid 2D function tabulation')
            n += 1
          # Left
          for q in range(numQuadPts-1, 0, -1):
            self.quadPts[n][0] = quadpts[0]
            self.quadPts[n][1] = quadpts[q]
            self.quadWts[n]    = quadwts[0]*quadwts[q]
            m = 0
            # Bottom
            for g in range(0, numBasisFns-1):
              self.basis[n][m] = basis[0][g]*basis[q][0]
              self.basisDeriv[n][m][0] = basisDeriv[0][g][0]*basis[q][0]
              self.basisDeriv[n][m][1] = basis[0][g]*basisDeriv[q][0][0]
              m += 1
            # Right
            for f in range(0, numBasisFns-1):
              self.basis[n][m] = basis[0][numBasisFns-1]*basis[q][f]
              self.basisDeriv[n][m][0] = basisDeriv[0][numBasisFns-1][0]*basis[q][f]
              self.basisDeriv[n][m][1] = basis[0][numBasisFns-1]*basisDeriv[q][f][0]
              m += 1
            # Top
            for g in range(numBasisFns-1, 0, -1):
              self.basis[n][m] = basis[0][g]*basis[q][numBasisFns-1]
              self.basisDeriv[n][m][0] = basisDeriv[0][g][0]*basis[q][numBasisFns-1]
              self.basisDeriv[n][m][1] = basis[0][g]*basisDeriv[q][numBasisFns-1][0]
              m += 1
            # Left
            for f in range(numBasisFns-1, 0, -1):
              self.basis[n][m] = basis[0][0]*basis[q][f]
              self.basisDeriv[n][m][0] = basisDeriv[0][0][0]*basis[q][f]
              self.basisDeriv[n][m][1] = basis[0][0]*basisDeriv[q][f][0]
              m += 1
            # Interior
            for f in range(1, numBasisFns-1):
              for g in range(1, numBasisFns-1):
                self.basis[n][m] = basis[0][g]*basis[0][f]
                self.basisDeriv[q][0][f][g][0] = basisDeriv[0][g][0]*basis[q][f]
                self.basisDeriv[q][0][f][g][1] = basis[0][g]*basisDeriv[q][f][0]
                m += 1
            if not m == self.numCorners: raise RuntimeError('Invalid 2D function tabulation')
            n += 1
          # Interior
          for q in range(1, numQuadPts-1):
            for r in range(1, numQuadPts-1):
              self.quadPts[n][0] = quadpts[r]
              self.quadPts[n][1] = quadpts[q]
              self.quadWts[n]    = quadwts[r]*quadwts[q]
              m = 0
              # Bottom
              for g in range(0, numBasisFns-1):
                self.basis[n][m] = basis[r][g]*basis[q][0]
                self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[q][0]
                self.basisDeriv[n][m][1] = basis[r][g]*basisDeriv[q][0][0]
                m += 1
              # Right
              for f in range(0, numBasisFns-1):
                self.basis[n][m] = basis[r][numBasisFns-1]*basis[q][f]
                self.basisDeriv[n][m][0] = basisDeriv[r][numBasisFns-1][0]*basis[q][f]
                self.basisDeriv[n][m][1] = basis[r][numBasisFns-1]*basisDeriv[q][f][0]
                m += 1
              # Top
              for g in range(numBasisFns-1, 0, -1):
                self.basis[n][m] = basis[r][g]*basis[q][numBasisFns-1]
                self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[q][numBasisFns-1]
                self.basisDeriv[n][m][1] = basis[r][g]*basisDeriv[q][numBasisFns-1][0]
                m += 1
              # Left
              for f in range(numBasisFns-1, 0, -1):
                self.basis[n][m] = basis[r][0]*basis[q][f]
                self.basisDeriv[n][m][0] = basisDeriv[r][0][0]*basis[q][f]
                self.basisDeriv[n][m][1] = basis[r][0]*basisDeriv[q][f][0]
                m += 1
              # Interior
              for f in range(1, numBasisFns-1):
                for g in range(1, numBasisFns-1):
                  self.basis[n][m] = basis[r][g]*basis[q][f]
                  self.basisDeriv[q][r][f][g][0] = basisDeriv[r][g][0]*basis[q][f]
                  self.basisDeriv[q][r][f][g][1] = basis[r][g]*basisDeriv[q][f][0]
                  m += 1
              if not m == self.numCorners: raise RuntimeError('Invalid 2D function tabulation')
              n += 1
          if not n == self.numQuadPts: raise RuntimeError('Invalid 2D quadrature')
        elif dim == 3:
          self.vertices = numpy.zeros((self.numCorners, dim))
          n = 0
          # Depth
          for s in range(numBasisFns):
            # Bottom
            for r in range(0, numBasisFns-1):
              self.vertices[n][0] = vertices[r]
              self.vertices[n][1] = vertices[0]
              self.vertices[n][2] = vertices[s]
              n += 1
            # Right
            for q in range(0, numBasisFns-1):
              self.vertices[n][0] = vertices[numBasisFns-1]
              self.vertices[n][1] = vertices[q]
              self.vertices[n][2] = vertices[s]
              n += 1
            # Top
            for r in range(numBasisFns-1, 0, -1):
              self.vertices[n][0] = vertices[r]
              self.vertices[n][1] = vertices[numBasisFns-1]
              self.vertices[n][2] = vertices[s]
              n += 1
            # Left
            for q in range(numBasisFns-1, 0, -1):
              self.vertices[n][0] = vertices[0]
              self.vertices[n][1] = vertices[q]
              self.vertices[n][2] = vertices[s]
              n += 1
            # Interior
            for q in range(1, numBasisFns-1):
              for r in range(1, numBasisFns-1):
                self.vertices[n][0] = vertices[r]
                self.vertices[n][1] = vertices[q]
                self.vertices[n][2] = vertices[s]
                n += 1
          if not n == self.numCorners: raise RuntimeError('Invalid 3D function tabulation: '+str(n)+' should be '+str(self.numCorners))

          self.quadPts    = numpy.zeros((numQuadPts*numQuadPts*numQuadPts, dim))
          self.quadWts    = numpy.zeros((numQuadPts*numQuadPts*numQuadPts,))
          self.basis      = numpy.zeros((numQuadPts*numQuadPts*numQuadPts,
                                         numBasisFns*numBasisFns*numBasisFns))
          self.basisDeriv = numpy.zeros((numQuadPts*numQuadPts*numQuadPts,
                                         numBasisFns*numBasisFns*numBasisFns,
                                         dim))
          n = 0
          # Depth
          for s in range(numQuadPts):
            # Bottom
            for r in range(0, numQuadPts-1):
              self.quadPts[n][0] = quadpts[r]
              self.quadPts[n][1] = quadpts[0]
              self.quadPts[n][2] = quadpts[s]
              self.quadWts[n]    = quadwts[r]*quadwts[0]*quadwts[s]
              m = 0
              for h in range(numBasisFns):
                # Bottom
                for g in range(0, numBasisFns-1):
                  self.basis[n][m] = basis[r][g]*basis[0][0]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[0][0]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[r][g]*basisDeriv[0][0][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[r][g]*basis[0][0]*basisDeriv[s][h][0]
                  m += 1
                # Right
                for f in range(0, numBasisFns-1):
                  self.basis[n][m] = basis[r][numBasisFns-1]*basis[0][f]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[r][numBasisFns-1][0]*basis[0][f]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[r][numBasisFns-1]*basisDeriv[0][f][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[r][numBasisFns-1]*basis[0][f]*basisDeriv[s][h][0]
                  m += 1
                # Top
                for g in range(numBasisFns-1, 0, -1):
                  self.basis[n][m] = basis[r][g]*basis[0][numBasisFns-1]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[0][numBasisFns-1]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[r][g]*basisDeriv[0][numBasisFns-1][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[r][g]*basis[0][numBasisFns-1]*basisDeriv[s][h][0]
                  m += 1
                # Left
                for f in range(numBasisFns-1, 0, -1):
                  self.basis[n][m] = basis[r][0]*basis[0][f]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[r][0][0]*basis[0][f]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[r][0]*basisDeriv[0][f][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[r][0]*basis[0][f]*basisDeriv[s][h][0]
                  m += 1
                # Interior
                for f in range(1, numBasisFns-1):
                  for g in range(1, numBasisFns-1):
                    self.basis[n][m] = basis[r][g]*basis[0][f]*basis[s][h]
                    self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[0][f]*basis[s][h]
                    self.basisDeriv[m][m][1] = basis[r][g]*basisDeriv[0][f][0]*basis[s][h]
                    self.basisDeriv[n][m][2] = basis[r][g]*basis[0][f]*basisDeriv[s][h][0]
                    m += 1
              if not m == self.numCorners: raise RuntimeError('Invalid 3D function tabulation')
              n += 1
            # Right
            for q in range(0, numQuadPts-1):
              self.quadPts[n][0] = quadpts[numQuadPts-1]
              self.quadPts[n][1] = quadpts[q]
              self.quadPts[n][2] = quadpts[s]
              self.quadWts[n]    = quadwts[numQuadPts-1]*quadwts[q]*quadwts[s]
              m = 0
              for h in range(numBasisFns):
                # Bottom
                for g in range(0, numBasisFns-1):
                  self.basis[n][m] = basis[numQuadPts-1][g]*basis[q][0]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[numQuadPts-1][g][0]*basis[q][0]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[numQuadPts-1][g]*basisDeriv[q][0][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[numQuadPts-1][g]*basis[q][0]*basisDeriv[s][h][0]
                  m += 1
                # Right
                for f in range(0, numBasisFns-1):
                  self.basis[n][m] = basis[numQuadPts-1][numBasisFns-1]*basis[q][f]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[numQuadPts-1][numBasisFns-1][0]*basis[q][f]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[numQuadPts-1][numBasisFns-1]*basisDeriv[q][f][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[numQuadPts-1][numBasisFns-1]*basis[q][f]*basisDeriv[s][h][0]
                  m += 1
                # Top
                for g in range(numBasisFns-1, 0, -1):
                  self.basis[n][m] = basis[numQuadPts-1][g]*basis[q][numBasisFns-1]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[numQuadPts-1][g][0]*basis[q][numBasisFns-1]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[numQuadPts-1][g]*basisDeriv[q][numBasisFns-1][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[numQuadPts-1][g]*basis[q][numBasisFns-1]*basisDeriv[s][h][0]
                  m += 1
                # Left
                for f in range(numBasisFns-1, 0, -1):
                  self.basis[n][m] = basis[numQuadPts-1][0]*basis[q][f]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[numQuadPts-1][0][0]*basis[q][f]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[numQuadPts-1][0]*basisDeriv[q][f][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[numQuadPts-1][0]*basis[q][f]*basisDeriv[s][h][0]
                  m += 1
                # Interior
                for f in range(1, numBasisFns-1):
                  for g in range(1, numBasisFns-1):
                    self.basis[n][m] = basis[numQuadPts-1][g]*basis[q][f]*basis[s][h]
                    self.basisDeriv[n][m][0] = basisDeriv[numQuadPts-1][g][0]*basis[q][f]*basis[s][h]
                    self.basisDeriv[m][m][1] = basis[numQuadPts-1][g]*basisDeriv[q][f][0]*basis[s][h]
                    self.basisDeriv[n][m][2] = basis[numQuadPts-1][g]*basis[q][f]*basisDeriv[s][h][0]
                    m += 1
              if not m == self.numCorners: raise RuntimeError('Invalid 3D function tabulation')
              n += 1
            # Top
            for r in range(numQuadPts-1, 0, -1):
              self.quadPts[n][0] = quadpts[r]
              self.quadPts[n][1] = quadpts[numQuadPts-1]
              self.quadPts[n][2] = quadpts[s]
              self.quadWts[n]    = quadwts[r]*quadwts[numQuadPts-1]*quadwts[s]
              m = 0
              for h in range(numBasisFns):
                # Bottom
                for g in range(0, numBasisFns-1):
                  self.basis[n][m] = basis[r][g]*basis[numQuadPts-1][0]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[numQuadPts-1][0]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[r][g]*basisDeriv[numQuadPts-1][0][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[r][g]*basis[numQuadPts-1][0]*basisDeriv[s][h][0]
                  m += 1
                # Right
                for f in range(0, numBasisFns-1):
                  self.basis[n][m] = basis[r][numBasisFns-1]*basis[numQuadPts-1][f]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[r][numBasisFns-1][0]*basis[numQuadPts-1][f]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[r][numBasisFns-1]*basisDeriv[numQuadPts-1][f][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[r][numBasisFns-1]*basis[numQuadPts-1][f]*basisDeriv[s][h][0]
                  m += 1
                # Top
                for g in range(numBasisFns-1, 0, -1):
                  self.basis[n][m] = basis[r][g]*basis[numQuadPts-1][numBasisFns-1]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[numQuadPts-1][numBasisFns-1]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[r][g]*basisDeriv[numQuadPts-1][numBasisFns-1][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[r][g]*basis[numQuadPts-1][numBasisFns-1]*basisDeriv[s][h][0]
                  m += 1
                # Left
                for f in range(numBasisFns-1, 0, -1):
                  self.basis[n][m] = basis[r][0]*basis[numQuadPts-1][f]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[r][0][0]*basis[numQuadPts-1][f]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[r][0]*basisDeriv[numQuadPts-1][f][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[r][0]*basis[numQuadPts-1][f]*basisDeriv[s][h][0]
                  m += 1
                # Interior
                for f in range(1, numBasisFns-1):
                  for g in range(1, numBasisFns-1):
                    self.basis[n][m] = basis[r][g]*basis[numQuadPts-1][f]*basis[s][h]
                    self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[numQuadPts-1][f]*basis[s][h]
                    self.basisDeriv[m][m][1] = basis[r][g]*basisDeriv[numQuadPts-1][f][0]*basis[s][h]
                    self.basisDeriv[n][m][2] = basis[r][g]*basis[numQuadPts-1][f]*basisDeriv[s][h][0]
                    m += 1
              if not m == self.numCorners: raise RuntimeError('Invalid 3D function tabulation')
              n += 1
            # Left
            for q in range(numQuadPts-1, 0, -1):
              self.quadPts[n][0] = quadpts[0]
              self.quadPts[n][1] = quadpts[q]
              self.quadPts[n][2] = quadpts[s]
              self.quadWts[n]    = quadwts[0]*quadwts[q]*quadwts[s]
              m = 0
              for h in range(numBasisFns):
                # Bottom
                for g in range(0, numBasisFns-1):
                  self.basis[n][m] = basis[0][g]*basis[q][0]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[0][g][0]*basis[q][0]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[0][g]*basisDeriv[q][0][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[0][g]*basis[q][0]*basisDeriv[s][h][0]
                  m += 1
                # Right
                for f in range(0, numBasisFns-1):
                  self.basis[n][m] = basis[0][numBasisFns-1]*basis[q][f]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[0][numBasisFns-1][0]*basis[q][f]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[0][numBasisFns-1]*basisDeriv[q][f][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[0][numBasisFns-1]*basis[q][f]*basisDeriv[s][h][0]
                  m += 1
                # Top
                for g in range(numBasisFns-1, 0, -1):
                  self.basis[n][m] = basis[0][g]*basis[q][numBasisFns-1]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[0][g][0]*basis[q][numBasisFns-1]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[0][g]*basisDeriv[q][numBasisFns-1][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[0][g]*basis[q][numBasisFns-1]*basisDeriv[s][h][0]
                  m += 1
                # Left
                for f in range(numBasisFns-1, 0, -1):
                  self.basis[n][m] = basis[0][0]*basis[q][f]*basis[s][h]
                  self.basisDeriv[n][m][0] = basisDeriv[0][0][0]*basis[q][f]*basis[s][h]
                  self.basisDeriv[n][m][1] = basis[0][0]*basisDeriv[q][f][0]*basis[s][h]
                  self.basisDeriv[n][m][2] = basis[0][0]*basis[q][f]*basisDeriv[s][h][0]
                  m += 1
                # Interior
                for f in range(1, numBasisFns-1):
                  for g in range(1, numBasisFns-1):
                    self.basis[n][m] = basis[0][g]*basis[q][f]*basis[s][h]
                    self.basisDeriv[n][m][0] = basisDeriv[0][g][0]*basis[q][f]*basis[s][h]
                    self.basisDeriv[m][m][1] = basis[0][g]*basisDeriv[q][f][0]*basis[s][h]
                    self.basisDeriv[n][m][2] = basis[0][g]*basis[q][f]*basisDeriv[s][h][0]
                    m += 1
              if not m == self.numCorners: raise RuntimeError('Invalid 3D function tabulation')
              n += 1
            # Interior
            for q in range(1, numQuadPts-1):
              for r in range(1, numQuadPts-1):
                self.quadPts[n][0] = quadpts[r]
                self.quadPts[n][1] = quadpts[q]
                self.quadPts[n][2] = quadpts[s]
                self.quadWts[n]    = quadwts[r]*quadwts[q]*quadwts[s]
                m = 0
                for h in range(numBasisFns):
                  # Bottom
                  for g in range(0, numBasisFns-1):
                    self.basis[n][m] = basis[r][g]*basis[q][0]*basis[s][h]
                    self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[q][0]*basis[s][h]
                    self.basisDeriv[n][m][1] = basis[r][g]*basisDeriv[q][0][0]*basis[s][h]
                    self.basisDeriv[n][m][2] = basis[r][g]*basis[q][0]*basisDeriv[s][h][0]
                    m += 1
                  # Right
                  for f in range(0, numBasisFns-1):
                    self.basis[n][m] = basis[r][numBasisFns-1]*basis[q][f]*basis[s][h]
                    self.basisDeriv[n][m][0] = basisDeriv[r][numBasisFns-1][0]*basis[q][f]*basis[s][h]
                    self.basisDeriv[n][m][1] = basis[r][numBasisFns-1]*basisDeriv[q][f][0]*basis[s][h]
                    self.basisDeriv[n][m][2] = basis[r][numBasisFns-1]*basis[q][f]*basisDeriv[s][h][0]
                    m += 1
                  # Top
                  for g in range(numBasisFns-1, 0, -1):
                    self.basis[n][m] = basis[r][g]*basis[q][numBasisFns-1]*basis[s][h]
                    self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[q][numBasisFns-1]*basis[s][h]
                    self.basisDeriv[n][m][1] = basis[r][g]*basisDeriv[q][numBasisFns-1][0]*basis[s][h]
                    self.basisDeriv[n][m][2] = basis[r][g]*basis[q][numBasisFns-1]*basisDeriv[s][h][0]
                    m += 1
                  # Left
                  for f in range(numBasisFns-1, 0, -1):
                    self.basis[n][m] = basis[r][0]*basis[q][f]*basis[s][h]
                    self.basisDeriv[n][m][0] = basisDeriv[r][0][0]*basis[q][f]*basis[s][h]
                    self.basisDeriv[n][m][1] = basis[r][0]*basisDeriv[q][f][0]*basis[s][h]
                    self.basisDeriv[n][m][2] = basis[r][0]*basis[q][f]*basisDeriv[s][h][0]
                    m += 1
                  # Interior
                  for f in range(1, numBasisFns-1):
                    for g in range(1, numBasisFns-1):
                      self.basis[n][m] = basis[r][g]*basis[q][f]*basis[s][h]
                      self.basisDeriv[n][m][0] = basisDeriv[r][g][0]*basis[q][f]*basis[s][h]
                      self.basisDeriv[m][m][1] = basis[r][g]*basisDeriv[q][f][0]*basis[s][h]
                      self.basisDeriv[n][m][2] = basis[r][g]*basis[q][f]*basisDeriv[s][h][0]
                      m += 1
                if not m == self.numCorners: raise RuntimeError('Invalid 3D function tabulation')
          if not n == self.numQuadPts: raise RuntimeError('Invalid 2D quadrature')
        self.vertices = numpy.reshape(self.vertices, (self.numCorners, dim))
        self.quadPts = numpy.reshape(self.quadPts, (self.numQuadPts, dim))
        self.quadWts = numpy.reshape(self.quadWts, (self.numQuadPts))
        self.basis = numpy.reshape(self.basis, (self.numQuadPts, self.numCorners))
        self.basisDeriv = numpy.reshape(self.basisDeriv, (self.numQuadPts, self.numCorners, dim))
    else:
      # Need 0-D quadrature for boundary conditions of 1-D meshes
      self.cellDim = 0
      self.numCorners = 1
      self.numQuadPts = 1
      self.basis = numpy.array([1.0])
      self.basisDeriv = numpy.array([1.0])
      self.quadPts = numpy.array([0.0])
      self.quadWts = numpy.array([1.0])

    self._info.line("Cell geometry: ")
    self._info.line(self.geometry)
    self._info.line("Vertices: ")
    self._info.line(self.vertices)
    self._info.line("Quad pts:")
    self._info.line(quadrature.get_points())
    self._info.line("Quad wts:")
    self._info.line(quadrature.get_weights())
    self._info.line("Basis fns @ quad pts ):")
    self._info.line(self.basis)
    self._info.line("Basis fn derivatives @ quad pts:")
    self._info.line(self.basisDeriv)
    self._info.log()    
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    import FIAT.shapes

    ReferenceCell._configure(self)
    self.cellDim = self.inventory.dimension
    self.degree = self.inventory.degree
    self.order = self.inventory.order

    if self.order == -1:
      self.order = 2*self.degree+1
    return


  def _setupGeometry(self, spaceDim):
    """
    Setup reference cell geometry object.
    """
    self.geometry = None
    if 3 == self.cellDim:
      if 3 == spaceDim:
        from geometry.GeometryHex3D import GeometryHex3D
        self.geometry = GeometryHex3D()
    elif 2 == self.cellDim:
      if 2 == spaceDim:
        from geometry.GeometryQuad2D import GeometryQuad2D
        self.geometry = GeometryQuad2D()
      elif 3 == spaceDim:
        from geometry.GeometryQuad3D import GeometryQuad3D
        self.geometry = GeometryQuad3D()
    elif 1 == self.cellDim:
      if 1 == spaceDim:
        from geometry.GeometryLine1D import GeometryLine1D
        self.geometry = GeometryLine1D()
      elif 2 == spaceDim:
        from geometry.GeometryLine2D import GeometryLine2D
        self.geometry = GeometryLine2D()
      elif 3 == spaceDim:
        from geometry.GeometryLine3D import GeometryLine3D
        self.geometry = GeometryLine3D()
    elif 0 == self.cellDim:
      if 1 == spaceDim:
        from geometry.GeometryPoint1D import GeometryPoint1D
        self.geometry = GeometryPoint1D()
      elif 2 == spaceDim:
        from geometry.GeometryPoint2D import GeometryPoint2D
        self.geometry = GeometryPoint2D()
      elif 3 == spaceDim:
        from geometry.GeometryPoint3D import GeometryPoint3D
        self.geometry = GeometryPoint3D()
    if None == self.geometry:
      raise ValueError("Could not set shape of cell for '%s' in spatial " \
                       "dimension '%s'." % (self.name, spaceDim))
    return
  

  def _setupQuadrature(self):
    """
    Setup quadrature rule for reference cell.
    """
    import FIAT.quadrature
    import FIAT.shapes
    return FIAT.quadrature.make_quadrature_by_degree(FIAT.shapes.LINE, self.order)


  def _setupElement(self):
    """
    Setup the finite element for reference cell.
    """
    import FIAT.Lagrange
    import FIAT.shapes
    return FIAT.Lagrange.Lagrange(FIAT.shapes.LINE, self.degree)


  def _setupVertices(self, element):
    """
    Setup evaluation functional points for reference cell.
    """
    return element.Udual.pts


# FACTORIES ////////////////////////////////////////////////////////////

def reference_cell():
  """
  Factory associated with FIATLagrange.
  """
  return FIATLagrange()


# End of file 
