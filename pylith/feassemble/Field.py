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

## @file pyre/feassemble/Field.py
## @brief Python PyLith field.

from pyre.components.Component import Component

def validateFamilyOrder(value):
  raise NotImplementedError, "validateFamilyOrder() not implemented."
  return value


def validateQuadratureOrder(value):
  raise NotImplementedError, "validateQuadratureOrder() not implemented."
  return value


def validateShape(value):
  raise NotImplementedError, "validateShape() not implemented."
  return value


# Field class
class Field(Component):
  """Python finite-element assembler."""

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing Field facilities and properties."""

    ## @class Inventory
    ## Python object for managing Field facilities and properties.
    ##
    ## \b Properties
    ## @li \b family Element family for field
    ## @li \b family_order Order of element family
    ## @li \b quadrature_order Order for quadrature
    ## @li \b shape Element shape
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    family = pyre.inventory.str("family", default="lagrange")
    family.meta['tip'] = "Element family for field"

    familyOrder = pyre.inventory.int("family_order", default=1,
                                     validator=validateFamilyOrder)
    familyOrder.meta['tip'] = "Order of element family"

    quadratureOrder = pyre.inventory.int("quadrature_order", default=1,
                                         validator=validateQuadratureOrder)
    quadratureOrder.meta['tip'] = "Order for quadrature."

    shape = pyre.inventory.str("shape", default="tetrahedron",
                               validator=validateShape)
    shape.meta['tip'] = "Element shape."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def initialize(self):
    """Setup basis fns and quadrature info."""
    return


  def __init__(self, name="field"):
    """Constructor."""
    Component.__init__(self, name, facility="field")
    self.sieveField = None
    self.basisInfo = None
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """Set members based using inventory."""
    self.family = self.inventory.family
    self.familyOrder = self.inventory.familyOrder
    self.quadratureOrder = self.inventory.quadratureOrder
    self.shape = self.inventory.shape
    return
  

# version
__id__ = "$Id$"

# End of file 
