# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/topology/MeshImporter.py
#
# @brief Python implementation of importing a mesh.
#
# Factory: mesh_generator.

from .MeshGenerator import MeshGenerator


class MeshImporter(MeshGenerator):
    """Python implementation of importing a mesh.

    FACTORY: mesh_generator.
    """

    import pythia.pyre.inventory

    reorderMesh = pythia.pyre.inventory.bool("reorder_mesh", default=True)
    reorderMesh.meta['tip'] = "Reorder mesh using reverse Cuthill-McKee."

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    reader = pythia.pyre.inventory.facility("reader", family="mesh_io", factory=MeshIOAscii)
    reader.meta['tip'] = "Mesh reader."

    from .Distributor import Distributor
    distributor = pythia.pyre.inventory.facility("distributor", family="mesh_distributor", factory=Distributor)
    distributor.meta['tip'] = "Mesh distributor."

    from .MeshRefiner import MeshRefiner
    refiner = pythia.pyre.inventory.facility("refiner", family="mesh_refiner", factory=MeshRefiner)
    refiner.meta['tip'] = "Mesh refiner."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="meshimporter"):
        """Constructor.
        """
        MeshGenerator.__init__(self, name)
        self._loggingPrefix = "MeIm "
        return

    def preinitialize(self, problem):
        """Do minimal initialization.
        """
        MeshGenerator.preinitialize(self, problem)

        self.reader.preinitialize()
        self.distributor.preinitialize()
        self.refiner.preinitialize()
        return

    def create(self, problem, faults=None):
        """Hook for creating mesh.
        """
        from pylith.utils.profiling import resourceUsageString
        from pylith.mpi.Communicator import petsc_comm_world
        comm = petsc_comm_world()

        self._setupLogging()
        logEvent = "%screate" % self._loggingPrefix
        self._eventLogger.eventBegin(logEvent)

        # Read mesh
        mesh = self.reader.read(self.debug)
        if self.debug:
            mesh.view()

        # Reorder mesh
        if self.reorderMesh:
            logEvent2 = "%sreorder" % self._loggingPrefix
            self._eventLogger.eventBegin(logEvent2)
            self._debug.log(resourceUsageString())
            if 0 == comm.rank:
                self._info.log("Reordering cells and vertices.")
            from pylith.topology.ReverseCuthillMcKee import ReverseCuthillMcKee
            ordering = ReverseCuthillMcKee()
            ordering.reorder(mesh)
            self._eventLogger.eventEnd(logEvent2)

        # Adjust topology
        self._debug.log(resourceUsageString())
        if 0 == comm.rank:
            self._info.log("Adjusting topology.")
        self._adjustTopology(mesh, faults, problem)

        # Distribute mesh
        if comm.size > 1:
            if 0 == comm.rank:
                self._info.log("Distributing mesh.")
            self.distributor.initialize()
            mesh = self.distributor.distribute(mesh, problem.normalizer)
            if self.debug:
                mesh.view()
            mesh.memLoggingStage = "DistributedMesh"

        # Refine mesh (if necessary)
        newMesh = self.refiner.refine(mesh)
        if not newMesh == mesh:
            mesh.cleanup()
            newMesh.memLoggingStage = "RefinedMesh"

        # Can't reorder mesh again, because we do not have routine to
        # unmix normal and hybrid cells.

        # Nondimensionalize mesh (coordinates of vertices).
        from pylith.topology.topology import MeshOps_nondimensionalize
        MeshOps_nondimensionalize(newMesh, problem.normalizer)

        self._eventLogger.eventEnd(logEvent)
        return newMesh

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based on inventory.
        """
        MeshGenerator._configure(self)
        return

    def _setupLogging(self):
        """Setup event logging.
        """
        MeshGenerator._setupLogging(self)
        self._eventLogger.registerEvent("%sreorder" % self._loggingPrefix)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_generator():
    """Factory associated with MeshImporter.
    """
    return MeshImporter()


# End of file
