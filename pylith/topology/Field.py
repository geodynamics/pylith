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

from .topology import Field as ModuleField


class Field(ModuleField):
    """
    Vector field over the simulation domain or a piece of the domain.
    """

    def __init__(self, mesh):
        """Constructor.
        """
        ModuleField.__init__(self, mesh)

    def cleanup(self):
        """Deallocate PETSc and local data structures.
        """
        self.deallocate()


# End of file
