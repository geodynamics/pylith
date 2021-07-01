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
# @file pylith/topology/ReverseCuthillMcKee.py
#
# @brief Python interface to PETSc reverse Cuthill-McKee reordering
# of mesh cells and vertices.

from .topology import ReverseCuthillMcKee as ModuleReverseCuthillMcKee


class ReverseCuthillMcKee(ModuleReverseCuthillMcKee):
    """Python interface to PETSc reverse Cuthill-McKee reordering of mesh
    cells and vertices.
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self):
        """Constructor.
        """
        return

    def reorder(self, mesh):
        """Set communicator.
        """
        ModuleReverseCuthillMcKee.reorder(mesh)
        return


# End of file
