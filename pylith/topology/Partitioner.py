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

## @file pylith/topology/Partitioner.py
##
## @brief Python manager for mesh partitioner.
##
## Factory: mesh_partitioner.

from pyre.components.Component import Component

# Partitioner class
class Partitioner(Component):
  """
  Python manager for mesh partitioner.

  Factory: mesh_partitioner
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Partitioner facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Partitioner facilities and properties.
    ##
    ## \b Properties
    ## @li \b partitioner Name of mesh partitioner {"parmetis", "chaco"}
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    partitioner = pyre.inventory.str("partitioner", default="chaco",
                                     validator=pyre.inventory.choice(["chaco",
                                                                      "parmetis"]))
    partitioner.meta['tip'] = "Name of mesh partitioner."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="partitioner"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="partitioner")
    import pylith.topology.topology as bindings
    self.cppHandle = bindings.Partitioner()
    return


  def distribute(self, mesh):
    """
    Distribute a Mesh
    """
    return self.cppHandle.distribute(mesh, self.partitioner)


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.partitioner = self.inventory.partitioner
    return


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_partitioner():
  """
  Factory associated with Partitioner.
  """
  return Partitioner()


# End of file 
