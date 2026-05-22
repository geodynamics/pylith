# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from pylith.utils.PetscComponent import PetscComponent


class InitializePhase(PetscComponent):
    """
    Abstract base class for mesh initialization phases.
    """

    def __init__(self, name="initialize_phase"):
        """Constructor."""
        PetscComponent.__init__(self, name, facility="initialize_phase")

    def preinitialize(self):
        """Preinitialize mesh initialization phase."""
        self._createModuleObj()

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)

    def _createModuleObj(self):
        raise NotImplementedError("Implement _createModuleObj().")

# End of file
