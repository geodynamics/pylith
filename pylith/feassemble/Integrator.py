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

## @file pylith/feassemble/Integrator.py

## @brief Python abstract base class for integration of operator
## actions with finite-elements.

from pyre.components.Component import Component

# Integrator class
class Integrator(Component):
  """
  Python abstract base class for integration of actions with
  finite-elements.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Integrator facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Integrator facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b quadrature Quadrature object for integration
    ## @li \b db Database for material properties.

    import pyre.inventory

    from Quadrature import Quadrature
    quadrature = pyre.inventory.facility("quadrature", factory=Quadrature)
    quadrature.meta['tip'] = "Quadrature object for integration."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="integrator"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="integrator")
    self.cppHandle = None
    self.quadrature = None
    return


  def initialize(self, mesh):
    """
    Initialize C++ integrator object.
    """
    q = self.quadrature
    q.initialize()
    self.cppHandle.quadrature = q.cppHandle
    self.cppHandle.createParameters(mesh.cppHandle)
    return
  
  
  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.quadrature = self.inventory.quadrature
    self.db = self.inventory.db
    return


# End of file 
