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
from .initializers import MeshReader as ModuleMeshReader


class MeshReader(InitializePhase, ModuleMeshReader):
    """
    Read mesh.

    Implements `InitializePhase`.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.mesh_initializer.phases.read_mesh]
            reader = pylith.meshio.MeshIOPetsc
        """
    }

    import pythia.pyre.inventory
    from pylith.meshio.MeshIOPetsc import MeshIOPetsc

    reader = pythia.pyre.inventory.facility(
        "reader", family="mesh_input", factory=MeshIOPetsc)
    reader.meta["tip"] = "Mesh reader."


    # PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////

    def __init__(self, name="mesh_reader"):
        """Constructor.
        """
        InitializePhase.__init__(self, name)

    def preinitialize(self):
        """Preinitialize.
        """
        InitializePhase.preinitialize(self)
        self.reader.preinitialize()
        ModuleMeshReader.setReader(self, self.reader)

    def _configure(self):
        """Set members based using inventory.
        """
        InitializePhase._configure(self)

    def _createModuleObj(self):
        """Create handle to C++ object."""
        ModuleMeshReader.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////
def initialize_phase():
    """Factory associated with MeshReader."""
    return MeshReader()


# End of file
