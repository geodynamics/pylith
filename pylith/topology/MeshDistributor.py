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

## @file pylith/topology/MeshDistributor.py
##
## @brief Python abstract base class for mesh distributor.
##
## Factory: mesh_distributor.

from pyre.components.Component import Component

# MeshDistributor class
class MeshDistributor(Component):
  """
  Python abstract base class for mesh distributor.

  Factory: mesh_distributor
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing MeshDistributor facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing MeshDistributor facilities and properties.
    ##
    ## \b Properties
    ## @li \b debug Debugging flag for mesh.
    ## @li \b interpolate Build intermediate mesh topology elements (if true)
    ## @li \b partitioner Name of mesh partitioner {"parmetis", "chaco"}
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    debug = pyre.inventory.bool("debug", default=False)
    debug.meta['tip'] = "Debugging flag for mesh."

    interpolate = pyre.inventory.bool("interpolate", default=False)
    interpolate.meta['tip'] = "Build intermediate mesh topology elements"

    partitioner = pyre.inventory.str("partitioner", default="chaco",
                                      validator=pyre.inventory.choice(["chaco",
                                                                       "parmetis"]))
    partitioner.meta['tip'] = "Name of mesh partitioner."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshdistributor"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh_distributor")
    import pylith.topology.topology as bindings
    self.cppHandle = bindings.MeshDistributor()
    self.debug = False
    self.interpolate = False
    self.partitioner = ''
    return


  def distributeMesh(self, mesh):
    """
    Distribute a Mesh
    """
    return self.cppHandle.distributeMesh(mesh, self.partitioner)


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.debug = self.inventory.debug
    self.interpolate = self.inventory.interpolate
    self.partitioner = self.inventory.partitioner
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
