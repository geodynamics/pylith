# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

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
