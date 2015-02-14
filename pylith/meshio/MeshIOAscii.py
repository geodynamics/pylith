#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pyre/meshio/MeshIOAscii.py
##
## @brief Python object for reading/writing finite-element mesh from
## simple ASCII file.
##
## Factory: mesh_io

from MeshIOObj import MeshIOObj
from meshio import MeshIOAscii as ModuleMeshIOAscii

# Validator for filename
def validateFilename(value):
  """
  Validate filename.
  """
  if 0 == len(value):
    msg = "Filename for ASCII input mesh not specified.\n\n" + \
        "PyLith runs some minimum validation on the input before even\n" + \
        "displaying help information. To test PyLith, run an example as\n" +\
        "discussed in the manual.\n"
    raise ValueError(msg)
  return value


# MeshIOAscii class
class MeshIOAscii(MeshIOObj, ModuleMeshIOAscii):
  """
  Python object for reading/writing finite-element mesh from simple
  ASCII file.

  Factory: mesh_io
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(MeshIOObj.Inventory):
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
    ## @li coordsys Coordinate system associated with mesh.

    import pyre.inventory

    filename = pyre.inventory.str("filename", default="", 
                                  validator=validateFilename)
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
    Create C++ MeshIOAscii object.
    """
    ModuleMeshIOAscii.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_io():
  """
  Factory associated with MeshIOAscii.
  """
  return MeshIOAscii()


# End of file 
