#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/materials/IsotropicPowerLaw.py
##
# @brief Python material for isotropic power-law viscoelastic material.
##
# Factory: elasticity_rheology

from .RheologyElasticity import RheologyElasticity
from .materials import IsotropicPowerLaw as ModuleLinearElasticity


class IsotropicPowerLaw(RheologyElasticity, ModuleLinearElasticity):
    """Python material for isotropic power-law viscoelastic material.

    FACTORY: elasticity_rheology
    """

    import pythia.pyre.inventory

    useReferenceState = pythia.pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="isotropicpowerlaw"):
        """Constructor.
        """
        RheologyElasticity.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsIsotropicPowerLaw import AuxSubfieldsIsotropicPowerLaw
        self.auxiliarySubfields = AuxSubfieldsIsotropicPowerLaw("auxiliary_subfields")

    def preinitialize(self, problem):
        RheologyElasticity.preinitialize(self, problem)

        ModuleLinearElasticity.useReferenceState(self, self.useReferenceState)

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleLinearElasticity.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def elasticity_rheology():
    """Factory associated with IsotropicPowerLaw.
    """
    return IsotropicPowerLaw()


# End of file
