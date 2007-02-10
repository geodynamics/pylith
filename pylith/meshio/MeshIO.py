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

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing MeshIO facilities and properties."""

    ## @class Inventory
    ## Python object for managing MeshIO facilities and properties.
    ##
    ## \b Properties
    ## @li \b interpolate Build intermediate mesh topology elements (if true)
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
    self.interpolate = False
    return


  def read(self):
    """
    Read finite-element mesh and store in Sieve mesh object.

    @returns PETSc mesh object containing finite-element mesh
    """
    from pylith.topology.Mesh import Mesh
    print "Creating Mesh object"
    mesh = Mesh()
    print "Setting interpolate"
    self.cppHandle.interpolate = self.interpolate
    print "Reading mesh"
    mesh.cppHandle = self.cppHandle.read(mesh.cppHandle)
    return 


  def write(self, mesh):
    """
    Write finite-element mesh.stored in Sieve mesh object.

    @param mesh PETSc mesh object containing finite-element mesh
    """
    self.cppHandle.interpolate = self.interpolate
    self.cppHandle.write(mesh.handle)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.interpolate = self.inventory.interpolate
    return


# End of file 
