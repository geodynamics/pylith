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

## @file pyre/meshio/MeshIOAscii.py
##
## @brief Python object for reading/writing finite-element mesh from
## simple ASCII file.
##
## Factory: mesh_io

from MeshIO import MeshIO

# MeshIOAscii class
class MeshIOAscii(MeshIO):
  """
  Python object for reading/writing finite-element mesh from simple
  ASCII file.

  Factory: mesh_io
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(MeshIO.Inventory):
    """
    Python object for managing MeshIOAscii facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing MeshIOAscii facilities and properties.
    ##
    ## \b Properties
    ## @li \b filename Name of mesh file
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    filename = pyre.inventory.str("filename", default="")
    filename.meta['tip'] = "Name of mesh file"

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pyre.inventory.facility("coordsys", family="coordsys",
                                       factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshioascii"):
    """
    Constructor.
    """
    MeshIO.__init__(self, name)
    import pylith.meshio.meshio as bindings
    self.cppHandle = bindings.MeshIOAscii()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    MeshIO._configure(self)
    self.filename = self.inventory.filename
    self.coordsys = self.inventory.coordsys
    return


  def _sync(self):
    """
    Force synchronization between Python and C++.
    """
    MeshIO._sync(self)
    self.cppHandle.filename = self.filename
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_io():
  """
  Factory associated with MeshIOAscii.
  """
  return MeshIOAscii()


# End of file 
