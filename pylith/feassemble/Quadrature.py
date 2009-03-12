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
# QuadratureBase class
class QuadratureBase(Component):
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
    return


  def preinitialize(self, spaceDim):
    """
    Setup quadrature object.
    """
    self._info.log("Initializing reference cell.")
    cell = self.cell
    cell.initialize(spaceDim)

    self._info.log("Initializing C++ quadrature.")
    self._initialize(cell)
    self.refGeometry(cell.geometry)
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
    self.minJacobian(self.inventory.minJacobian)
    self.checkConditioning(self.inventory.checkConditioning)
    self.cell = self.inventory.cell
    return


# ----------------------------------------------------------------------
from feassemble import MeshQuadrature as ModuleMeshQuadrature

# MeshQuadrature class
class MeshQuadrature(QuadratureBase, ModuleMeshQuadrature):
  """
  Python object for integrating over finite-elements using quadrature.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshquadrature"):
    """
    Constructor.
    """
    QuadratureBase.__init__(self, name)
    ModuleMeshQuadrature.__init__(self)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _initialize(self, cell):
    """
    Initialize C++ quadrature object.
    """
    print "cell.basis shape: ", cell.basis.shape
    print "cell.basisDeriv shape: ", cell.basisDeriv.shape
    print "cell.quadPts shape: ", cell.quadPts.shape
    print "cell.quadWts shape: ", cell.quadWts.shape
    print "cell.geometry.spaceDim: ", cell.geometry.spaceDim()
    ModuleMeshQuadrature.initialize(self, cell.basis, cell.basisDeriv,
                                    cell.quadPts, cell.quadWts,
                                    cell.geometry.spaceDim())
    return


# End of file 
