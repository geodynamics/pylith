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

## @file pylith/topology/MeshGenerator.py
##
## @brief Python abstract base class for mesh generator.
##
## Factory: mesh_generator.

from pyre.components.Component import Component

# MeshGenerator class
class MeshGenerator(Component):
  """
  Python abstract base class for mesh generator.

  Factory: mesh_generator
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing MeshGenerator facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing MeshGenerator facilities and properties.
    ##
    ## \b Properties
    ## @li \b debug Debugging flag for mesh.
    ## @li \b interpolate Build intermediate mesh topology elements (if true)
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    debug = pyre.inventory.bool("debug", default=False)
    debug.meta['tip'] = "Debugging flag for mesh."

    interpolate = pyre.inventory.bool("interpolate", default=False)
    interpolate.meta['tip'] = "Build intermediate mesh topology elements"


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshgenerator"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh_generator")
    self.debug = False
    self.interpolate = False
    return


  def create(self, faults=None):
    """
    Generate a Mesh from a boundary
    """
    raise NotImplementedError("MeshGenerator.create() not implemented.")
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.debug = self.inventory.debug
    self.interpolate = self.inventory.interpolate
    return


  def _adjustTopology(self, mesh, faults):
    """
    Adjust topology for fault implementation.
    """
    if not faults is None:
      for fault in faults:
        mesh.adjustTopology(fault)
    return
  

# End of file 
