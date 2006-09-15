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

## @brief Python object for reading/writing finite-element mesh from
## simple ASCII file.

from MeshIO import MeshIO

# MeshIOAscii class
class MeshIOAscii(MeshIO):
  """
  Python object for reading/writing finite-element mesh from simple
  ASCII file.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing MeshIOAscii facilities and properties."""

    ## @class Inventory
    ## Python object for managing MeshIOAscii facilities and properties.
    ##
    ## \b Properties
    ## @li \b filename Name of mesh file
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    filename = pyre.inventory.str("filename", default=False)
    filename.meta['tip'] = "Name of mesh file"


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshioascii"):
    """Constructor."""
    MeshIO.__init__(self, name)
    import pylith.meshio.meshio as bindings
    self.cppHandler = bindings.MeshIOAscii()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """Set members based using inventory."""
    MeshIO.configure(self)
    self.filename = self.inventory.filename
    return


# version
__id__ = "$Id$"

# End of file 
