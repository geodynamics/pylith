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

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshio"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh_io")
    self.cppHandle = None
    self.coordsys = None
    return


  def read(self, debug, interpolate):
    """
    Read finite-element mesh and store in Sieve mesh object.

    @returns PETSc mesh object containing finite-element mesh
    """
    self._info.log("Reading finite-element mesh")

    # Set flags
    self._sync()
    self.cppHandle.debug = debug
    self.cppHandle.interpolate = interpolate

    # Initialize coordinate system
    if self.coordsys is None:
      raise ValueError, "Coordinate system for mesh is unknown."
    self.coordsys.initialize()

    from pylith.topology.Mesh import Mesh
    mesh = Mesh()
    mesh.initialize(self.coordsys)

    # Read mesh
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
    return


  def _sync(self):
    """
    Force synchronization between Python and C++.
    """
    return


# End of file 
