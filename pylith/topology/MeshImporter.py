# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .MeshGenerator import MeshGenerator


class MeshImporter(MeshGenerator):
    """
    Base class for reading a finite-element mesh from files.

    Implements `MeshGenerator`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.meshimporter]
            reorder_mesh = True
            check_topology = True
            reader = pylith.meshio.MeshIOCubit
            refiner = pylith.topology.RefineUniform
        """
    }

    import pythia.pyre.inventory

    reorderMesh = pythia.pyre.inventory.bool("reorder_mesh", default=True)
    reorderMesh.meta['tip'] = "Reorder mesh using reverse Cuthill-McKee."

    checkTopology = pythia.pyre.inventory.bool("check_topology", default=True)
    checkTopology.meta['tip'] = "Check topology of imported mesh."

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    reader = pythia.pyre.inventory.facility("reader", family="mesh_io", factory=MeshIOAscii)
    reader.meta['tip'] = "Reader for mesh files."

    from .Distributor import Distributor
    distributor = pythia.pyre.inventory.facility("distributor", family="mesh_distributor", factory=Distributor)
    distributor.meta['tip'] = "Distributes mesh among processes."

    from .MeshRefiner import MeshRefiner
    refiner = pythia.pyre.inventory.facility("refiner", family="mesh_refiner", factory=MeshRefiner)
    refiner.meta['tip'] = "Performs uniform global mesh refinement after distribution among processes (default is no refinement)."

    def __init__(self, name="meshimporter"):
        """Constructor.
        """
        MeshGenerator.__init__(self, name)
        self._loggingPrefix = "PL.MeshImporter."

    def preinitialize(self, problem):
        """Do minimal initialization.
        """
        MeshGenerator.preinitialize(self, problem)

        self.reader.preinitialize()
        self.distributor.preinitialize()
        self.refiner.preinitialize()

    def create(self, problem, faults=None):
        """Hook for creating mesh.
        """
        from pylith.utils.profiling import resourceUsageString
        from pylith.mpi.Communicator import mpi_is_root
        isRoot = mpi_is_root()

        self._setupLogging()
        logEvent = f"{self._loggingPrefix}create"
        self._eventLogger.eventBegin(logEvent)

        # Read mesh
        mesh = self.reader.read(self.checkTopology)

        # Reorder mesh
        if self.reorderMesh:
            logEvent2 = f"{self._loggingPrefix}reorder"
            self._eventLogger.eventBegin(logEvent2)
            self._debug.log(resourceUsageString())
            if isRoot:
                self._info.log("Reordering cells and vertices.")
            from pylith.topology.ReverseCuthillMcKee import ReverseCuthillMcKee
            ordering = ReverseCuthillMcKee()
            ordering.reorder(mesh)
            self._eventLogger.eventEnd(logEvent2)

        # Adjust topology
        self._debug.log(resourceUsageString())
        if isRoot:
            self._info.log("Adjusting topology.")
        self._adjustTopology(mesh, faults, problem)

        # Distribute mesh
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if comm.size > 1:
            if isRoot:
                self._info.log("Distributing mesh.")
            mesh = self.distributor.distribute(mesh, problem)
            mesh.memLoggingStage = "DistributedMesh"

        # Refine mesh (if necessary)
        newMesh = self.refiner.refine(mesh)
        if not newMesh == mesh:
            mesh.cleanup()

        # Can't reorder mesh again, because we do not have routine to
        # unmix normal and hybrid cells.

        # Nondimensionalize mesh (coordinates of vertices).
        from pylith.topology.topology import MeshOps_nondimensionalize
        MeshOps_nondimensionalize(newMesh, problem.normalizer)

        self._eventLogger.eventEnd(logEvent)
        return newMesh

    def _configure(self):
        """Set members based on inventory.
        """
        MeshGenerator._configure(self)

    def _setupLogging(self):
        """Setup event logging.
        """
        MeshGenerator._setupLogging(self)
        self._eventLogger.registerEvent(f"{self._loggingPrefix}reorder")


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_generator():
    """Factory associated with MeshImporter.
    """
    return MeshImporter()


# End of file
