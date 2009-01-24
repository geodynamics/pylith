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

from MeshIOObj import MeshIOObj
from meshio import MeshIOCubit as ModuleMeshIOCubit

# MeshIOCubit class
class MeshIOCubit(MeshIOObj, ModuleMeshIOCubit):
  """
  Python object for reading/writing finite-element mesh from Cubit.

  Factory: mesh_io
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(MeshIOObj.Inventory):
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
    MeshIOObj.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    MeshIOObj._configure(self)
    self.coordsys = self.inventory.coordsys
    self.filename(self.inventory.filename)
    return


  def _createModuleObj(self):
    """
    Create C++ MeshIOCubit object.
    """
    ModuleMeshIOCubit.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_io():
  """
  Factory associated with MeshIOCubit.
  """
  return MeshIOCubit()


# End of file 
