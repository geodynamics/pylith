# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
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
            use_edge_weighting = True
            write_partition = False
            data_writer = pylith.meshio.DataWriterHDF5
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
        ModuleDistributor.setPartitioner(self, self.partitioner)
        ModuleDistributor.setUseEdgeWeighting(self, self.useEdgeWeighting)
        if self.writePartition:
            ModuleDistributor.setDataWriter(self, self.dataWriter)

    def _configure(self):
        """Set members using inventory.
        """
        PetscComponent._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_distributor():
    """Factory associated with Distributor.
    """
    return Distributor()


# End of file
