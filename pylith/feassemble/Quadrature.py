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

## @file pylith/feassemble/quadrature/Qudrature.py
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
    ## @li \b check_conditoning Check element matrices for 
    ##   ill-conditioning.
    ##
    ## \b Facilities
    ## @li \b cell Reference cell with basis functions and quadrature rules

    import pyre.inventory

    minJacobian = pyre.inventory.float("min_jacobian", default=1.0e-06)
    minJacobian.meta['tip'] = "Minimum allowable determinant of Jacobian."

    checkConditioning = pyre.inventory.bool("check_conditioning",
                                            default=False)
    checkConditioning.meta['tip'] = \
        "Check element matrices for ill-conditioning."

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


  def preinitialize(self):
    """
    Setup quadrature object.
    """
    self._createCppHandle()
    
    self.cppHandle.minJacobian = self.minJacobian
    self.cppHandle.checkConditioning = self.checkConditioning

    self._info.log("Initializing reference cell.")
    cell = self.cell
    cell.initialize(self.spaceDim)

    if cell.cellDim != self.cellDim:
      raise TypeError("Dimension of reference cell '%d' does not match "
                      "dimension of quadrature implementation '%d'." % \
                      (cell.cellDim, self.cellDim))


    self._info.log("Initializing C++ quadrature.")
    self.cppHandle.initialize(cell.basis, cell.basisDeriv,
                              cell.quadPts, cell.quadWts,
                              cell.cellDim, cell.numCorners, cell.numQuadPts,
                              self.spaceDim)
    self.cppHandle.refGeometry = cell.geometry.cppHandle
    return


  def initialize(self):
    """
    Initialize quadrature object.
    """
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.minJacobian = self.inventory.minJacobian
    self.checkConditioning = self.inventory.checkConditioning
    self.cell = self.inventory.cell
    return


  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    raise NotImplementedError("Please implement _createCppHandle() in " \
                              "derived class.")
  
  
# End of file 
