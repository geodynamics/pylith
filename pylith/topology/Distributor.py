# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.utils.PetscComponent import PetscComponent
from .topology import Distributor as ModuleDistributor


class Distributor(PetscComponent, ModuleDistributor):
    """
    Distributor of the the mesh among processes.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.mesh_generator.distributor]
            partitioner = parmetis
        """
    }

    import pythia.pyre.inventory

    partitioner = pythia.pyre.inventory.str("partitioner", default="parmetis",
                                     validator=pythia.pyre.inventory.choice(["parmetis", "chaco", "simple"]))
    partitioner.meta['tip'] = "Name of mesh partitioner (PETSc must be built with partitioner)."

    useEdgeWeighting = pythia.pyre.inventory.bool("use_edge_weighting", default=True)
    useEdgeWeighting.meta["tip"] = "Use edge weighting (parmetis only)."

    writePartition = pythia.pyre.inventory.bool("write_partition", default=False)
    writePartition.meta['tip'] = "Write partition information to file."

    from pylith.meshio.DataWriterHDF5 import DataWriterHDF5
    dataWriter = pythia.pyre.inventory.facility("data_writer", factory=DataWriterHDF5, family="data_writer")
    dataWriter.meta['tip'] = "Data writer for partition information."

    def __init__(self, name="mesh_distributor"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="mesh_distributor")

    def preinitialize(self):
        """Do minimal initialization."""
        ModuleDistributor.__init__(self)

    def distribute(self, mesh, problem):
        """Distribute a Mesh
        """
        self._setupLogging()
        logEvent = f"{self._loggingPrefix}distribute"
        self._eventLogger.eventBegin(logEvent)

        from pylith.topology.Mesh import Mesh
        newMesh = Mesh(mesh.getDimension())
        ModuleDistributor.distribute(newMesh, mesh, problem.interfaces.components(), self.partitioner, self.useEdgeWeighting)

        mesh.cleanup()

        if self.writePartition:
            self.dataWriter.setFilename(problem.defaults.outputDir, problem.defaults.simName, "partition")
            ModuleDistributor.write(self.dataWriter, newMesh)

        self._eventLogger.eventEnd(logEvent)
        return newMesh

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)

    def _setupLogging(self):
        """Setup event logging.
        """
        self._loggingPrefix = "PL.Distributor."
        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.setClassName("Distributor")
        logger.initialize()
        events = ["distribute"]
        for event in events:
            logger.registerEvent(f"{self._loggingPrefix}{event}")

        self._eventLogger = logger


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_distributor():
    """Factory associated with Distributor.
    """
    return Distributor()


# End of file
