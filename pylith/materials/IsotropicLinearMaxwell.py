# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .RheologyElasticity import RheologyElasticity
from .materials import IsotropicLinearMaxwell as ModuleLinearElasticity


class IsotropicLinearMaxwell(RheologyElasticity, ModuleLinearElasticity):
    """
    Isotropic linear Maxwell viscoelastic bulk rheology.

    Implements `RheologyElasticity`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_maxwell.rheology]
            use_reference_state = False

            auxiliary_subfields.shear_modulus.basis_order = 0
            auxiliary_subfields.bulk_modulus.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    useReferenceState = pythia.pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    def __init__(self, name="isotropiclinearmaxwell"):
        """Constructor.
        """
        RheologyElasticity.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsIsotropicLinearMaxwell import AuxSubfieldsIsotropicLinearMaxwell
        self.auxiliarySubfields = AuxSubfieldsIsotropicLinearMaxwell("auxiliary_subfields")

    def preinitialize(self, problem):
        RheologyElasticity.preinitialize(self, problem)

        ModuleLinearElasticity.useReferenceState(self, self.useReferenceState)

        return

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleLinearElasticity.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def elasticity_rheology():
    """Factory associated with IsotropicLinearMaxwell.
    """
    return IsotropicLinearMaxwell()


# End of file
