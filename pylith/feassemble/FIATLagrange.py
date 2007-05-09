b#!/usr/bin/env python
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
## Factory: reference_cell.

from FIATCell import FIATCell

import numpy

def validateDimension(self, dim):
  if dim < 1 or dim > 3:
    raise ValueError("Dimension of Lagrange element must be 1, 2, or 3.")
  return

# FIATLagrange class
class FIATLagrange(FIATCell):
  """
  Python object for managing basis functions and quadrature rules of a
  Lagrange reference finite-element cell using FIAT.

  Factory: reference_cell.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(FIATCell.Inventory):
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

    quadOrder = pyre.inventory.int("quad_order", default=-1)
    quadOrder.meta['tip'] = "Order of quadrature rule."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fiatlagrange"):
    """
    Constructor.
    """
    FIATCell.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    import FIAT.shapes

    FIATCell._configure(self)
    self.cellDim = self.inventory.dimension
    self.degree = self.inventory.degree
    self.quadOrder = self.inventory.quadOrder
    if self.cellDim == 1:
      self.shape = FIAT.shapes.LINE
    elif self.cellDim == 2:
      self.shape = FIAT.shapes.TRIANGLE
    elif self.cellDim == 3:
      self.shape = FIAT.shapes.TETRAHEDRON
    if self.quadOrder == -1:
      self.quadOrder = 2*self.degree+1
    self.numCorners = self.cellDim+1
    return


  def _setupQuadrature(self):
    """
    Setup quadrature rule for reference cell.
    """
    import FIAT.quadrature
    self.quadrature = FIAT.quadrature.make_quadrature_by_degree(shape, self.quadOrder)
    self.numQuadPts = len(quadrature.get_points())
    self.quadPts = quadrature.get_points()
    self.quadWts = quadrature.get_weights()
    return


  def _setupBasisFns(self):
    """
    Setup basis functions for reference cell.
    """
    from FIAT.Lagrange import Lagrange

    self.element = Lagrange(self.shape, self.degree)
    points = self.quadrature.get_points()
    basis = self.element.function_space()
    self.numBasisFuncs = len(basis)
    self.basis = numpy.transpose(basis.tabulate(points))
    self.basisDeriv = numpy.transpose([basis.deriv_all(d).tabulate(points) for d in range(self.cellDim)])
    return

# FACTORIES ////////////////////////////////////////////////////////////

def reference_cell():
  """
  Factory associated with FIATLagrange.
  """
  return FIATLagrange()


# End of file 
