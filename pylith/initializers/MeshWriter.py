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
from .initializers import MeshWriter as ModuleMeshWriter


class MeshWriter(InitializePhase):
    """
    Write mesh.

    Implements `InitializePhase`.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.mesh_initializer.phases.write_mesh]
            writer = pylith.meshio.MeshIOPetsc
        """
    }

    import pythia.pyre.inventory
    from pylith.meshio.MeshIOPetsc import MeshIOPetsc

    writer = pythia.pyre.inventory.facility(
        "writer", family="mesh_output", factory=MeshIOPetsc)
    writer.meta["tip"] = "Mesh writer."


    # PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////

    def __init__(self, name="mesh_writer"):
        """Constructor.
        """
        InitializePhase.__init__(self, name)

    def preinitialize(self):
        """Preinitialize.
        """
        InitializePhase.preinitialize(self)
        self.writer.preinitialize()
        ModuleMeshWriter.setWriter(self, self.writer)

    def _configure(self):
        """Set members based using inventory.
        """
        InitializePhase._configure(self)

    def _createModuleObj(self):
        """Create handle to C++ object."""
        ModuleMeshWriter.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////
def initialize_phase():
    """Factory associated with MeshWriter."""
    return MeshWriter()


# End of file
