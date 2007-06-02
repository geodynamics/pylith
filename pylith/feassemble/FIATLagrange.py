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

def validateDimension(self, dim):
  if dim < 1 or dim > 3:
    raise ValueError("Dimension of Lagrange element must be 1, 2, or 3.")
  return

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


  def initialize(self):
    """
    Initialize reference finite-element cell from a tensor product of
    1-D Lagrange elements.
    """
    quadrature = self._setupQuadrature()
    element    = self._setupElement()
    dim        = self.cellDim
    
    # Evaluate basis functions at quadrature points
    quadpts     = numpy.array(quadrature.get_points())
    quadwts     = numpy.array(quadrature.get_weights())
    numQuadPts  = len(quadpts)
    basis       = numpy.array(element.function_space().tabulate(quadrature.get_points())).transpose()
    numBasisFns = len(element.function_space())

    # Evaluate derivatives of basis functions at quadrature points
    basisDeriv = numpy.array([element.function_space().deriv_all(d).tabulate(quadrature.get_points()) \
                              for d in range(1)]).transpose()

    self.numQuadPts = numQuadPts**dim
    self.numCorners = numBasisFns**dim

    if dim == 1:
      self.quadPts    = quadpts
      self.quadWts    = quadwts
      self.basis      = basis
      self.basisDeriv = basisDeriv
    else:
      if dim == 2:
        self.quadPts    = numpy.zeros((numQuadPts, numQuadPts, dim))
        self.quadWts    = numpy.zeros((numQuadPts, numQuadPts))
        self.basis      = numpy.zeros((numQuadPts, numQuadPts, numBasisFns, numBasisFns))
        self.basisDeriv = numpy.zeros((numQuadPts, numQuadPts, numBasisFns, numBasisFns, dim))
        for q in range(numQuadPts):
          for r in range(numQuadPts):
            self.quadPts[q][r][0] = quadpts[q]
            self.quadPts[q][r][1] = quadpts[r]
            self.quadWts[q][r]    = quadwts[q]*quadwts[r]
            for f in range(numBasisFns):
              for g in range(numBasisFns):
                self.basis[q][r][f][g]         = basis[q][f]*basis[r][g]
                self.basisDeriv[q][r][f][g][0] = basisDeriv[q][f][0]*basis[r][g]
                self.basisDeriv[q][r][f][g][1] = basis[q][f]*basisDeriv[r][g][0]
      elif dim == 3:
        self.quadPts    = numpy.zeros((numQuadPts, numQuadPts, numQuadPts, dim))
        self.quadWts    = numpy.zeros((numQuadPts, numQuadPts, numQuadPts))
        self.basis      = numpy.zeros((numQuadPts, numQuadPts, numQuadPts, numBasisFns, numBasisFns, numBasisFns))
        self.basisDeriv = numpy.zeros((numQuadPts, numQuadPts, numQuadPts, numBasisFns, numBasisFns, numBasisFns, dim))
        for q in range(numQuadPts):
          for r in range(numQuadPts):
            for s in range(numQuadPts):
              self.quadPts[q][r][s][0] = quadpts[q]
              self.quadPts[q][r][s][1] = quadpts[r]
              self.quadPts[q][r][s][2] = quadpts[s]
              self.quadWts[q][r][s]    = quadwts[q]*quadwts[r]*quadwts[s]
              for f in range(numBasisFns):
                for g in range(numBasisFns):
                  for h in range(numBasisFns):
                    self.basis[q][r][s][f][g][h]         = basis[q][f]*basis[r][g]*basis[s][h]
                    self.basisDeriv[q][r][s][f][g][h][0] = basisDeriv[q][f][0]*basis[r][g]*basis[s][h]
                    self.basisDeriv[q][r][s][f][g][h][1] = basis[q][f]*basisDeriv[r][g][0]*basis[s][h]
                    self.basisDeriv[q][r][s][f][g][h][2] = basis[q][f]*basis[r][g]*basisDeriv[s][h][0]
      self.quadPts    = numpy.reshape(self.quadPts,    (self.numQuadPts, dim))
      self.quadWts    = numpy.reshape(self.quadWts,    (self.numQuadPts))
      self.basis      = numpy.reshape(self.basis,      (self.numQuadPts, self.numCorners))
      self.basisDeriv = numpy.reshape(self.basisDeriv, (self.numQuadPts, self.numCorners, dim))
    self._info.line("Basis (quad pts):")
    self._info.line(self.basis)
    self._info.line("Basis derivatives (quad pts):")
    self._info.line(self.basisDeriv)
    self._info.line("Quad pts:")
    self._info.line(self.quadPts)
    self._info.line("Quad wts:")
    self._info.line(self.quadWts)
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


# FACTORIES ////////////////////////////////////////////////////////////

def reference_cell():
  """
  Factory associated with FIATLagrange.
  """
  return FIATLagrange()


# End of file 
