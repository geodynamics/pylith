# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from pylith.petsc.Component import Component
from .meshio import OutputTrigger as ModuleOutputTrigger


class OutputTrigger(Component, ModuleOutputTrigger):
    """
    Abstract base class for managing how often output is written.
    """

    def __init__(self, name="outputtrigger"):
        """Constructor.
        """
        Component.__init__(self, name, facility="outputtrigger")

    def _configure(self):
        """Set members using inventory.
        """
        Component._configure(self)


# End of file
