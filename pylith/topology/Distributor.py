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
# @file pylith/topology/Distributor.py
#
# @brief Python manager for distributing mesh among processors.
#
# Factory: mesh_distributor.

from pylith.utils.PetscComponent import PetscComponent
from .topology import Distributor as ModuleDistributor


class Distributor(PetscComponent, ModuleDistributor):
    """Python manager for distributing mesh among processors.

    FACTORY: mesh_distributor
    """

    import pythia.pyre.inventory

    partitioner = pythia.pyre.inventory.str("partitioner", default="chaco",
                                     validator=pythia.pyre.inventory.choice(["chaco", "metis", "parmetis", "simple"]))
    partitioner.meta['tip'] = "Name of mesh partitioner."

    writePartition = pythia.pyre.inventory.bool("write_partition", default=False)
    writePartition.meta['tip'] = "Write partition information to file."

    from pylith.meshio.DataWriterVTK import DataWriterVTK
    dataWriter = pythia.pyre.inventory.facility("data_writer", factory=DataWriterVTK, family="data_writer")
    dataWriter.meta['tip'] = "Data writer for partition information."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="mesh_distributor"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="mesh_distributor")
        return

    def preinitialize(self):
        """Do minimal initialization."""
        ModuleDistributor.__init__(self)
        return

    def distribute(self, mesh, normalizer):
        """Distribute a Mesh
        """
        self._setupLogging()
        logEvent = "%sdistribute" % self._loggingPrefix
        self._eventLogger.eventBegin(logEvent)

        from pylith.topology.Mesh import Mesh
        newMesh = Mesh(mesh.dimension())
        if self.partitioner == "metis":
            partitionerName = "parmetis"
        else:
            partitionerName = self.partitioner
        ModuleDistributor.distribute(newMesh, mesh, partitionerName)

        mesh.cleanup()

        if self.writePartition:
            self.dataWriter.initialize(normalizer)
            ModuleDistributor.write(self.dataWriter, newMesh)

        self._eventLogger.eventEnd(logEvent)
        return newMesh

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)
        return

    def _setupLogging(self):
        """Setup event logging.
        """
        self._loggingPrefix = "Dist "
        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.setClassName("FE Distribution")
        logger.initialize()
        events = ["distribute"]
        for event in events:
            logger.registerEvent("%s%s" % (self._loggingPrefix, event))

        self._eventLogger = logger
        return


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_distributor():
    """Factory associated with Distributor.
    """
    return Distributor()


# End of file
