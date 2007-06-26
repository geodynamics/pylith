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

## @file pylith/topology/MeshImporter.py
##
## @brief Python implementation of importing a mesh.
##
## Factory: mesh_generator.

from MeshGenerator import MeshGenerator

# MeshImporter class
class MeshImporter(MeshGenerator):
  """
  Python implementation of importing a mesh.

  Factory: mesh_generator.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(MeshGenerator.Inventory):
    """
    Python object for managing MeshImporter facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing MeshImporter facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b importer Mesh importer.
    ## @li \b partitioner Mesh partitioner.

    import pyre.inventory

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = pyre.inventory.facility("importer", family="mesh_io",
                                       factory=MeshIOAscii)
    importer.meta['tip'] = "Mesh importer."

    from Distributor import Distributor
    distributor = pyre.inventory.facility("distributor",
                                          family="mesh_distributor",
                                          factory=Distributor)
    distributor.meta['tip'] = "Mesh distributor."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshimporter"):
    """
    Constructor.
    """
    MeshGenerator.__init__(self, name)
    return


  def create(self, faults=None):
    """
    Hook for creating mesh.
    """
    mesh = self.importer.read(self.debug, self.interpolate)
    self._info.log("Adjusting topology.")
    self._adjustTopology(mesh, faults)
    self._info.log("Distributing mesh.")
    mesh = self.distributor.distribute(mesh)
    mesh.view()
    return mesh


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based on inventory.
    """
    MeshGenerator._configure(self)
    self.importer = self.inventory.importer
    self.distributor = self.inventory.distributor
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_generator():
  """
  Factory associated with MeshImporter.
  """
  return MeshImporter()


# End of file 
