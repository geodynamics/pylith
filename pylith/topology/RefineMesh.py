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


class RefineMesh(PetscComponent):
    """
    Abstract base class for refining a mesh in parallel.
    """

    def __init__(self, name="refiner"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="refiner")

    def preinitialize(self):
        """Do minimal initialization."""
        self._createModuleObj()

    def _configure(self):
        """Set members using inventory.
        """
        PetscComponent._configure(self)

    def _createModuleObj(self):
        raise NotImplementedError("Implement _createModuleObj().")


# End of file
