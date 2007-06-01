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
    basisFns = self._setupBasisFns()
    
    # Evaluate basis functions at quadrature points
    quadpts = quadrature.get_points()
    basis = numpy.array(basisFns.tabulate(quadpts)).transpose()
    self.basis = numpy.reshape(basis.flatten(), basis.shape)

    # Evaluate derivatives of basis functions at quadrature points
    import FIAT.shapes
    dim = FIAT.shapes.dimension(basisFns.base.shape)
    basisDeriv = numpy.array([basisFns.deriv_all(d).tabulate(quadpts) \
                              for d in range(dim)]).transpose()
    self.basisDeriv = numpy.reshape(basisDeriv.flatten(), basisDeriv.shape)

    self.quadPts = numpy.array(quadrature.get_points())
    self.quadWts = numpy.array(quadrature.get_weights())

    self.numCorners = len(basisFns)
    self.numQuadPts = len(quadrature.get_weights())

    self._info.line("Basis (quad pts):")
    self._info.line(self.basis)
    self._info.line("Basis derivatives (quad pts):")
    self._info.line(self.basisDeriv)
    self._info.line("Quad pts:")
    self._info.line(quadrature.get_points())
    self._info.line("Quad wts:")
    self._info.line(quadrature.get_weights())
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
    # :TODO: ADD STUFF HERE
    return


  def _setupBasisFns(self):
    """
    Setup basis functions for reference cell.
    """
    from FIAT.Lagrange import Lagrange

    # :TODO: ADD STUFF HERE
    return


# FACTORIES ////////////////////////////////////////////////////////////

def reference_cell():
  """
  Factory associated with FIATLagrange.
  """
  return FIATLagrange()


# End of file 
