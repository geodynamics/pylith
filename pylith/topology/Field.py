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
# @file pylith/topology/Field.py
#
# @brief Python object for managing a vector field over vertices or
# cells of a finite-element mesh.

from .topology import Field as ModuleField


class Field(ModuleField):
    """Python object for managing a vector field over vertices or cells of
    a finite-element mesh.
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, mesh):
        """Constructor.
        """
        ModuleField.__init__(self, mesh)
        return

    def cleanup(self):
        """Deallocate PETSc and local data structures.
        """
        self.deallocate()
        return


# End of file
