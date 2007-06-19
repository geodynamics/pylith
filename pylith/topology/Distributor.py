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

## @file pylith/topology/Distributor.py
##
## @brief Python manager for distributing mesh among processors.
##
## Factory: mesh_distributor.

from pyre.components.Component import Component

# Distributor class
class Distributor(Component):
  """
  Python manager for distributing mesh among processors.

  Factory: mesh_distributor
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Distributor facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Distributor facilities and properties.
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
    self.cppHandle = None
    return


  def distribute(self, mesh):
    """
    Distribute a Mesh
    """
    self._createCppHandle()
    
    from Mesh import Mesh
    newMesh = Mesh()
    newMesh.cppHandle = self.cppHandle.distribute(mesh.cppHandle,
                                                  self.partitioner)
    newMesh.coordsys = mesh.coordsys
    return newMesh


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.partitioner = self.inventory.partitioner
    return


  def _createCppHandle(self):
    """
    Create handle to C++ object.
    """
    if None == self.cppHandle:
      import pylith.topology.topology as bindings
      self.cppHandle = bindings.Distributor()
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_partitioner():
  """
  Factory associated with Distributor.
  """
  return Distributor()


# End of file 
