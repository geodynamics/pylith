# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


from pylith.petsc.Component import Component


class EmptyBin(Component):
    """
    Empty container for a collection of objects.
    """

    def __init__(self, name="emptybin"):
        """Constructor."""
        Component.__init__(self, name, facility="empty_bin")


# End of file
