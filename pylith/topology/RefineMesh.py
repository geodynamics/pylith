# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


from pylith.petsc.Component import Component


class RefineMesh(Component):
    """
    Abstract base class for refining a mesh in parallel.
    """

    def __init__(self, name="refiner"):
        """Constructor."""
        Component.__init__(self, name, facility="refiner")

    def preinitialize(self):
        """Do minimal initialization."""
        self._createModuleObj()

    def _configure(self):
        """Set members using inventory."""
        Component._configure(self)

    def _createModuleObj(self):
        raise NotImplementedError("Implement _createModuleObj().")


# End of file
