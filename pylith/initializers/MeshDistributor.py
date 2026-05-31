# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .InitializePhase import InitializePhase
from .initializers import MeshDistributor as ModuleMeshDistributor


class MeshDistributor(InitializePhase):
    """
    Distribute mesh.

    Implements `InitializePhase`.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.mesh_initializer.phases.distribute_mesh]
            distributor = pylith.topology.Distributor
        """
    }

    import pythia.pyre.inventory
    from pylith.topology.Distributor import Distributor

    distributor = pythia.pyre.inventory.facility(
        "distributor", family="mesh_distributor", factory=Distributor)
    distributor.meta["tip"] = "Mesh distributor."


    # PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////

    def __init__(self, name="mesh_refiner"):
        """Constructor.
        """
        InitializePhase.__init__(self, name)

    def preinitialize(self):
        """Preinitialize.
        """
        InitializePhase.preinitialize(self)
        self.distributor.preinitialize()
        ModuleMeshDistributor.setDistributor(self, self.distributor)

    def _configure(self):
        """Set members using inventory.
        """
        InitializePhase._configure(self)

    def _createModuleObj(self):
        """Create handle to C++ object."""
        ModuleMeshDistributor.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////
def initialize_phase():
    """Factory associated with MeshDistributor."""
    return MeshDistributor()


# End of file
