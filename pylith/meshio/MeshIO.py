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
##
## @brief Python abstract base class for finite-element mesh I/O.
##
## Factory: mesh_io

from pyre.components.Component import Component

# MeshIO class
class MeshIO(Component):
  """
  Python abstract base class for finite-element mesh I/O.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing MeshIO facilities and properties.

    Factory: mesh.
    """

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
    """
    Component.__init__(self, name, facility="mesh_io")
    self.cppHandle = None
    self.interpolate = False
    self.coordsys = None
    return


  def read(self):
    """
    Read finite-element mesh and store in Sieve mesh object.

    @returns PETSc mesh object containing finite-element mesh
    """
    self._info.log("Reading finite-element mesh")
    self._sync()
    from pylith.topology.Mesh import Mesh
    mesh = Mesh()
    if self.coordsys is None:
      raise ValueError, "Coordinate system for mesh is unknown."
    mesh.initialize(self.coordsys)
    self.cppHandle.read(mesh.cppHandle)
    return mesh


  def write(self, mesh):
    """
    Write finite-element mesh.stored in Sieve mesh object.

    @param mesh PETSc mesh object containing finite-element mesh
    """
    self._info.log("Writing finite-element mesh")
    self._sync()
    self.cppHandle.write(mesh.cppHandle)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.interpolate = self.inventory.interpolate
    return


  def _sync(self):
    """
    Force synchronization between Python and C++.
    """
    self.cppHandle.interpolate = self.interpolate
    return
  

# End of file 
