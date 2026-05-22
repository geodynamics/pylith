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
from .initializers import MeshRefiner as ModuleMeshRefiner


class MeshRefiner(InitializePhase):
    """
    Refine mesh.

    Implements `InitializePhase`.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.mesh_initializer.phases.refine_mesh]
            refiner = pylith.topology.RefineUniform
        """
    }

    import pythia.pyre.inventory
    from pylith.topology.RefineUniform import RefineUniform

    refiner = pythia.pyre.inventory.facility(
        "refiner", family="refine", factory=RefineUniform)
    refiner.meta["tip"] = "Mesh refiner."


    # PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////

    def __init__(self, name="mesh_refiner"):
        """Constructor.
        """
        InitializePhase.__init__(self, name)

    def preinitialize(self):
        """Preinitialize.
        """
        InitializePhase.preinitialize(self)
        self.refiner.preinitialize()
        ModuleMeshRefiner.setRefiner(self, self.refiner)

    def _configure(self):
        """Set members based using inventory.
        """
        InitializePhase._configure(self)

    def _createModuleObj(self):
        """Create handle to C++ object."""
        ModuleMeshRefiner.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////
def initialize_phase():
    """Factory associated with MeshRefiner."""
    return MeshRefiner()


# End of file
