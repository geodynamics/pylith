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

## @file pylith/feassemble/Qudrature.py
##
## @brief Python abstract base class for integrating over
## finite-elements using quadrature.
##
## Factory: quadrature.

from pyre.components.Component import Component

# ----------------------------------------------------------------------
# Quadrature class
class Quadrature(Component):
  """
  Python abstract base class for integrating over finite-elements
  using quadrature.

  Factory: quadrature.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing Quadrature facilities and properties."""

    ## @class Inventory
    ## Python object for managing Quadrature facilities and properties.
    ##
    ## \b Properties
    ## @li \b min_jacobian Minimum allowable determinant of Jacobian.
    ##
    ## \b Facilities
    ## @li \b cell Reference cell with basis functions and quadrature rules

    import pyre.inventory

    minJacobian = pyre.inventory.float("min_jacobian", default=1.0e-06)
    minJacobian.meta['tip'] = "Minimum allowable determinant of Jacobian."

    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = pyre.inventory.facility("cell", family="reference_cell",
                                   factory=FIATSimplex)
    cell.meta['tip'] = "Reference cell with basis fns and quadrature rules."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="quadrature")
    self.minJacobian = 1.0e-06
    self.cppHandle = None
    self.spaceDim = None
    return


  def initialize(self):
    """
    Initialize C++ quadrature object.
    """
    if self.cppHandle is None:
      raise ValueError("C++ handle not set.")
    
    self.cppHandle.minJacobian = self.minJacobian

    self._info.log("Initializing reference cell.")
    c = self.cell
    c.initialize()

    if c.cellDim != self.cellDim:
      raise TypeError("Dimension of reference cell '%d' does not match "
                      "dimension of quadrature implementation '%d'." % \
                      (c.cellDim, self.cellDim))


    self._info.log("Initializing C++ quadrature.")
    self.cppHandle.initialize(c.basisVert, c.basisDerivVert,
                              c.basisQuad, c.basisDerivQuad,
                              c.quadPts, c.quadWts,
                              c.cellDim, c.numCorners, c.numQuadPts,
                              self.spaceDim)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.minJacobian = self.inventory.minJacobian
    self.cell = self.inventory.cell
    return


# End of file 
