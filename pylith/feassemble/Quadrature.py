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

## @file pyre/feassemble/Qudrature.py

## @brief Python abstract base class for integrating over
## finite-elements using quadrature.

# TODO
#
# Use FIAT to create numpy arrays containing basis functions and their
# derivatives evaludated at the quadrature points and the coordinates
# and weights of the quadrature points in the reference cell.

# DESIGN QUESTION
#
# The dimension of the space associated with the coordinates of
# the vertices is specific to the quadrature object. However, the
# quadrature points/weights and basis functions are associated with
# the reference cell not the actual cell.
#
# Where should the dimension of the space associated with the
# coordinates of the vertices be specified? I don't think the user
# will specify it, so that means it won't be in the inventory.

from pyre.components.Component import Component

# Quadrature class
class Quadrature(Component):
  """
  Python abstract base class for integrating over finite-elements
  using quadrature.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing Quadrature facilities and properties."""

    ## @class Inventory
    ## Python object for managing Quadrature facilities and properties.
    ##
    ## \b Properties
    ## @li \b jacobian_tolerance Minimum allowable determinant of Jacobian.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    jacobianTol = pyre.inventory.float("jacobian_tolerance", default=1.0e-06)
    jacobianTol.meta['tip'] = "Minimum allowable determinant of Jacobian."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="quadrature")
    #import pylith.feassemble.feassemble as bindings
    #self.cppHandle = bindings.Quadrature()
    self.spaceDim = 3
    return

  def initialize(self):
    """
    Initialize C++ quadrature object.
    """
    # Set minimum allowable determinant of Jacobian
    #self.cppHandle.jacobianTol = self.jacobianTol

    # Get basis functions, quadrature points
    #self.cppHandle.initialize(LOTS OF ARGS)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.jacobianTol = self.inventory.jacobianTol
    return


# version
__id__ = "$Id$"

# End of file 
