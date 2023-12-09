# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .topology import ReverseCuthillMcKee as ModuleReverseCuthillMcKee


class ReverseCuthillMcKee(ModuleReverseCuthillMcKee):
    """
    Interface to PETSc reverse Cuthill-McKee reordering of mesh cells and vertices.
    """

    def __init__(self):
        """Constructor.
        """
        return

    def reorder(self, mesh):
        """Set communicator.
        """
        ModuleReverseCuthillMcKee.reorder(mesh)


# End of file
