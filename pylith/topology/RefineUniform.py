# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.utils.PetscComponent import PetscComponent
from .topology import RefineUniform as ModuleRefineUniform


class RefineUniform(PetscComponent, ModuleRefineUniform):
    """
    Uniform global mesh refinement in parallel.

    Implements `MeshRefiner`.
    """
    DOC_CONFIG = {
        "cfg": """
            # Refine mesh twice to reduce size of cell edges by a factor of 4.
            [pylithapp.mesh_generator.refiner]
            levels = 0
        """
    }

    import pythia.pyre.inventory

    levels = pythia.pyre.inventory.int("levels", default=0, validator=pythia.pyre.inventory.greaterEqual(0))
    levels.meta['tip'] = "Number of refinement levels."

    def __init__(self, name="refineuniform"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="mesh_refiner")

    def preinitialize(self):
        """Do minimal initialization."""
        self._createModuleObj()
        ModuleRefineUniform.setNumLevels(self, self.levels)

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleRefineUniform.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_refiner():
    """Factory associated with RefineUniform.
    """
    return RefineUniform()


# End of file
