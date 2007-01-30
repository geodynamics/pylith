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
## @brief Python implementation of importing a mesh.

from MeshGenerator import MeshGenerator

# MeshImporter class
class MeshImporter(MeshGenerator):
  """
  Python implementation of importing a mesh.
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
    ## @li \b importer Mesh importer

    import pyre.inventory

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = pyre.inventory.facility("importer", factory=MeshIOAscii)
    importer.meta['tip'] = "Mesh importer."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def create(self):
    """
    Hook for creating mesh.
    """
    return self.importer.read()


  def __init__(self, name="meshimporter"):
    """
    Constructor.
    """
    MeshGenerator.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based on inventory.
    """
    self.importer = self.inventory.importer
    return
  

# End of file 
