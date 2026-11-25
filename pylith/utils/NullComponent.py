# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


from pylith.petsc.Component import Component


class NullComponent(Component):
    """
    Empty Pyre component.
    """

    def __init__(self):
        """Constructor."""
        Component.__init__(self, name="nullcomponent", facility="nullcomponent")

    def _cleanup(self):
        """Deallocate locally managed data structures."""
        return


# End of file
