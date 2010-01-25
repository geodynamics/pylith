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

from pylith.utils.PetscComponent import PetscComponent
from meshio import MeshIO as ModuleMeshIO

# MeshIOObj class
class MeshIOObj(PetscComponent, ModuleMeshIO):
  """
  Python abstract base class for finite-element mesh I/O.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshio"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="mesh_io")
    self.coordsys = None
    self._createModuleObj()
    return


  def read(self, debug, interpolate):
    """
    Read finite-element mesh and store in Sieve mesh object.

    @returns PETSc mesh object containing finite-element mesh
    """
    self._info.log("Reading finite-element mesh")

    # Set flags
    self.debug(debug)
    self.interpolate(interpolate)

    # Initialize coordinate system
    if self.coordsys is None:
      raise ValueError, "Coordinate system for mesh is unknown."

    from pylith.mpi.Communicator import petsc_comm_world
    from pylith.topology.Mesh import Mesh    
    mesh = Mesh()
    mesh.setComm(petsc_comm_world())
    mesh.coordsys(self.coordsys)

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
    PetscComponent._configure(self)
    return


  def _createModuleObj(self):
    """
    Create C++ MeshIO object.
    """
    raise NotImplementedError("MeshIO is an abstract base class.")
    return


# End of file 
