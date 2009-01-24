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

## @file pyre/meshio/MeshIOObj.py
##
## @brief Python abstract base class for finite-element mesh I/O.
##
## Factory: mesh_io

from pyre.components.Component import Component
from meshio import MeshIO as ModuleMeshIO

# MeshIOObj class
class MeshIOObj(Component, ModuleMeshIO):
  """
  Python abstract base class for finite-element mesh I/O.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshio"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh_io")
    self.coordsys = None
    self._createModuleObj()
    return


  def read(self, dim, normalizer, debug, interpolate):
    """
    Read finite-element mesh and store in Sieve mesh object.

    @returns PETSc mesh object containing finite-element mesh
    """
    self._info.log("Reading finite-element mesh")

    # Set flags
    self.normalizer(normalizer)
    self.debug(debug)
    self.interpolate(interpolate)

    # Initialize coordinate system
    if self.coordsys is None:
      raise ValueError, "Coordinate system for mesh is unknown."

    from pylith.mpi.Communicator import petsc_comm_world
    from pylith.topology.Mesh import Mesh    
    mesh = Mesh(petsc_comm_world(), dim)
    mesh.coordsys(self.coordsys)
    mesh.initialize()

    # Read mesh
    ModuleMeshIO.read(self, mesh)
    return mesh


  def write(self, mesh):
    """
    Write finite-element mesh.stored in Sieve mesh object.

    @param mesh PETSc mesh object containing finite-element mesh
    """
    self._info.log("Writing finite-element mesh")
    ModuleMeshIO.write(self, mesh)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    return


  def _createModuleObj(self):
    """
    Create C++ MeshIO object.
    """
    raise NotImplementedError("MeshIO is an abstract base class.")
    return


# End of file 
