# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

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
