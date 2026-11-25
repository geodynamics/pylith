# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


from pylith.petsc.Component import Component


class InitializePhase(Component):
    """
    Abstract base class for mesh initialization phases.
    """

    def __init__(self, name="initialize_phase"):
        """Constructor."""
        Component.__init__(self, name, facility="initialize_phase")

    def preinitialize(self):
        """Preinitialize mesh initialization phase."""
        self._createModuleObj()

    def _configure(self):
        """Set members using inventory."""
        Component._configure(self)

    def _createModuleObj(self):
        raise NotImplementedError("Implement _createModuleObj().")


# End of file
