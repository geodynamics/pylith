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
import FIAT.shapes

elementShapes = {'LINE':        FIAT.shapes.LINE,
                 'TRIANGLE':    FIAT.shapes.TRIANGLE,
                 'TETRAHEDRON': FIAT.shapes.TETRAHEDRON}

from pyre.components.Component import Component

def validateFamily(value):
  try:
    __import__('FIAT.'+str(value))
  except ImportError:
    raise ValueError, 'Invalid element family: '+str(value)
  return value

def validateShape(value):
  if not str(value).upper() in elementShapes:
    raise ValueError, 'Invalid element shape: '+str(value)
  return value.upper()

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

    family = pyre.inventory.str("family", default="Lagrange",
                                validator=validateFamily)
    family.meta['tip'] = "Element family for field"

    familyOrder = pyre.inventory.int("family_order", default=1,
                                     validator=pyre.inventory.greaterEqual(0))
    familyOrder.meta['tip'] = "Order of element family"

    quadratureOrder = pyre.inventory.int("quadrature_order", default=1,
                                         validator=pyre.inventory.greater(0))
    quadratureOrder.meta['tip'] = "Order for quadrature."

    shape = pyre.inventory.str("shape", default="tetrahedron",
                               validator=validateShape)
    shape.meta['tip'] = "Element shape."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def initialize(self):
    """Setup basis fns and quadrature info."""
    import FIAT.quadrature
    from pylith.utils import importing

    self._info.log('Creating the '+self.family+' element of order ' +
                   str(self.familyOrder))
    self.element = getattr(importing.importModule('FIAT.'+self.family),
                           self.family)(self.shape, self.familyOrder)
    self.quadrature = FIAT.quadrature.make_quadrature_by_degree(self.shape,
                                                        2*self.element.n - 1)
    self._info.log('Created quadrature of order '+str(self.quadrature.degree))
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
    self.shape = elementShapes[self.inventory.shape]
    return
  

# version
__id__ = "$Id$"

# End of file 
