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
from .initializers import MeshInsertInterfaces as ModuleMeshInsertInterfaces


class MeshInsertInterfaces(InitializePhase):
    """
    Insert fault interfaces using PETSc transform operation.

    Implements `InitializePhase`.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.mesh_initializer.phases.insert_interfaces]
            # No parameters
        """
    }

    def __init__(self, name="mesh_insert_interfaces"):
        """Constructor."""
        InitializePhase.__init__(self, name)

    def _configure(self):
        """Set members based on inventory."""
        InitializePhase._configure(self)

    def _createModuleObj(self):
        """Create handle to C++ object."""
        ModuleMeshInsertInterfaces.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////
def initialize_phase():
    """Factory associated with MeshInsertInterfaces."""
    return MeshInsertInterfaces()


# End of file
