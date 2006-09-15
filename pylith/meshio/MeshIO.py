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

## @file pyre/meshio/MeshIO.py

## @brief Python abstract base class for finite-element mesh I/O.

from pyre.components.Component import Component

# MeshIO class
class MeshIO(Component):
  """
  Python abstract base class for finite-element mesh I/O.
  """

  class Inventory(Component.Inventory):
    """Python object for managing MeshIO facilities and properties."""

    ## @class Inventory
    ## Python object for managing Field facilities and properties.
    ##
    ## \b Properties
    ## @li \b interpolate Build intermediate mesh topology elements
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    interpolate = pyre.inventory.bool("interpolate", default=False)
    interpolate.meta['tip'] = "Build intermediate mesh topology elements"


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshio"):
    """
    Constructor.

    @param name Component name
    """
    Component.__init__(self, name, facility="meshio")
    self.cppHandle = None
    return


  def read(self):
    """
    Read finite-element mesh and store in Sieve mesh object.

    @returns Sieve mesh object containing finite-element mesh
    """
    from pylith.topology.Mesh import Mesh
    mesh = Mesh()
    mesh.handle = self.cppHandle.read(self.interpolate)
    return 


  def write(self, mesh):
    """
    Write finite-element mesh.stored in Sieve mesh object.

    @param mesh Sieve mesh object containing finite-element mesh
    """
    self.cppHandle.write(mesh.handle)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """Set members based using inventory."""
    self.interpolate = self.inventory.interpolate
    return

# version
__id__ = "$Id$"

# End of file 
