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

## @brief Python object for managing basis functions and quadrature
## rules of a Lagrange reference finite-element cell using FIAT.

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

    quadOrder = pyre.inventory.int("quad_order", default=1)
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
    FIATCell._configure(self)
    self.dimension = self.inventory.dimension
    self.degree = self.inventory.degree
    self.order = self.inventory.order
    return


  def _setupQuadrature(self):
    """
    Setup quadrature rule for reference cell.
    """

    # ADD STUFF HERE
    raise NotImplementedError()

    return


  def _setupBasisFns(self):
    """
    Setup basis functions for reference cell.
    """
    from FIAT.Lagrange import Lagrange

    # ADD STUFF HERE
    raise NotImplementedError()

    return



# End of file 
