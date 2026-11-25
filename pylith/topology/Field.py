# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


from .topology import Field as ModuleField


class Field(ModuleField):
    """
    Vector field over the simulation domain or a piece of the domain.
    """

    def __init__(self, mesh):
        """Constructor."""
        ModuleField.__init__(self, mesh)

    def cleanup(self):
        """Deallocate PETSc and local data structures."""
        self.deallocate()


# End of file
