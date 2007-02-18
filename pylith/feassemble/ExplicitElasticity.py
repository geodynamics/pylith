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

## @file pylith/feassemble/ExplicitElasticity.py

## @brief Python object for explicit time integration of dynamic
## elasticity equation using finite-elements.

from IntegratorExplicit import IntegratorExplicit

# ExplicitElasticity class
class ExplicitElasticity(IntegratorExplicit):
  """
  Python object for explicit time integration of dynamic elasticity
  equation using finite-elements.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(IntegratorExplicit.Inventory):
    """
    Python object for managing ExplicitElasticity facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing ExplicitElasticity facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b db Database for material property parameters.

    import pyre.inventory

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = pyre.inventory.facility("db", factory=SimpleDB)
    db.meta['tip'] = "Database for material property parameters."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="explicitelasticity"):
    """
    Constructor.
    """
    IntegratorExplicit.__init__(self, name)

    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.ExplicitElasticity()
    return


  def initialize(self, mesh):
    """
    Initialize integrator.
    """
    self.cppHandle.setupMatProp(mesh.cppHandle, mesh.coordsys, db.cppHandle)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    IntegratorExplicit._configure(self)
    self.db = self.inventory.db
    return


# End of file 
