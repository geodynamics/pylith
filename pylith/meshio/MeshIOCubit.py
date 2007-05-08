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

## @file pyre/meshio/MeshIOCubit.py
##
## @brief Python object for reading/writing finite-element mesh from
## Cubit.
##
## Factory: mesh_io

from MeshIO import MeshIO

# MeshIOCubit class
class MeshIOCubit(MeshIO):
  """
  Python object for reading/writing finite-element mesh from Cubit.

  Factory: mesh_io
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(MeshIO.Inventory):
    """
    Python object for managing MeshIOCubit facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing MeshIOCubit facilities and properties.
    ##
    ## \b Properties
    ## @li \b filename Name of Cubit Exodus file.
    ##
    ## \b Facilities
    ## @li coordsys Coordinate system associated with mesh.

    import pyre.inventory

    filename = pyre.inventory.str("filename", default="mesh.exo")
    filename.meta['tip'] = "Name of Cubit Exodus file."

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pyre.inventory.facility("coordsys", family="coordsys",
                                       factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshiocubit"):
    """
    Constructor.
    """
    MeshIO.__init__(self, name)
    import pylith.meshio.meshio as bindings
    self.cppHandle = bindings.MeshIOCubit()
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
  Factory associated with MeshIOCubit.
  """
  return MeshIOCubit()


# End of file 
