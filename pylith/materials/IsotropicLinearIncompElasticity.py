# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .RheologyIncompressibleElasticity import RheologyIncompressibleElasticity
from .materials import IsotropicLinearIncompElasticity as ModuleLinearElasticity


class IsotropicLinearIncompElasticity(RheologyIncompressibleElasticity, ModuleLinearElasticity):
    """
    Isotropic linear incompressible elastic bulk rheology.

    Implements `RheologyIncompressibleElasticity`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_incompelastic.rheology]
            use_reference_state = False

            auxiliary_subfields.shear_modulus.basis_order = 0
            auxiliary_subfields.bulk_modulus.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    useReferenceState = pythia.pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    def __init__(self, name="isotropiclinearincompelasticity"):
        """Constructor.
        """
        RheologyIncompressibleElasticity.__init__(self, name)

    def _defaults(self):
        from .AuxSubfieldsIsotropicLinearElasticity import AuxSubfieldsIsotropicLinearElasticity
        self.auxiliarySubfields = AuxSubfieldsIsotropicLinearElasticity("auxiliary_subfields")

    def preinitialize(self, problem):
        RheologyIncompressibleElasticity.preinitialize(self, problem)
        ModuleLinearElasticity.useReferenceState(self, self.useReferenceState)

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleLinearElasticity.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def incompressible_elasticity_rheology():
    """Factory associated with IsotropicLinearIncompElasticity.
    """
    return IsotropicLinearIncompElasticity()


# End of file
