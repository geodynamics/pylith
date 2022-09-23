# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

from pylith.utils.PetscComponent import PetscComponent
from .meshio import MeshIO as ModuleMeshIO


class MeshIOObj(PetscComponent, ModuleMeshIO):
    """
    Abstract base class for finite-element mesh readers.
    """

    def __init__(self, name="meshio"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="mesh_io")
        self.coordsys = None

    def preinitialize(self):
        """Do minimal initialization."""
        self._createModuleObj()
        ModuleMeshIO.setIdentifier(self, self.aliases[-1])

    def read(self, debug):
        """Read finite-element mesh and store in Sieve mesh object.

        @returns PETSc DMPlex mesh object containing finite-element mesh
        """
        from pylith.mpi.Communicator import mpi_is_root
        if mpi_is_root():
            self._info.log("Reading finite-element mesh")

        # Set flags
        self.debug(debug)

        # Initialize coordinate system
        if self.coordsys is None:
            raise ValueError("Coordinate system for mesh is unknown.")

        from pylith.mpi.Communicator import petsc_comm_world
        from pylith.topology.Mesh import Mesh
        mesh = Mesh(dim=self.coordsys.getSpaceDim(), comm=petsc_comm_world())
        mesh.setCoordSys(self.coordsys)

        # Read mesh
        ModuleMeshIO.read(self, mesh)
        return mesh

    def write(self, mesh):
        """Write finite-element mesh.stored in Sieve mesh object.

        @param mesh PETSc mesh object containing finite-element mesh
        """
        from pylith.mpi.Communicator import mpi_is_root
        if mpi_is_root():
            self._info.log("Writing finite-element mesh")
        ModuleMeshIO.write(self, mesh)

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)

    def _createModuleObj(self):
        """Create C++ MeshIO object.
        """
        raise NotImplementedError("MeshIO is an abstract base class.")


# End of file
