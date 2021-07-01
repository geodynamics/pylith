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
# @file pylith/materials/IsotropicLinearIncompElasticity.py
#
# @brief Python material for isotropic, linearly elastic, incompressible
# material.
#
# Factory: incompressible_elasticity_rheology

from .RheologyIncompressibleElasticity import RheologyIncompressibleElasticity
from .materials import IsotropicLinearIncompElasticity as ModuleLinearElasticity


class IsotropicLinearIncompElasticity(RheologyIncompressibleElasticity, ModuleLinearElasticity):
    """Python material for isotropic, linearly elastic incompressible.

    FACTORY: incompressible_elasticity_rheology
    """

    import pythia.pyre.inventory

    useReferenceState = pythia.pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="isotropiclinearincompelasticity"):
        """Constructor.
        """
        RheologyIncompressibleElasticity.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsIsotropicLinearElasticity import AuxSubfieldsIsotropicLinearElasticity
        self.auxiliarySubfields = AuxSubfieldsIsotropicLinearElasticity("auxiliary_subfields")

    def preinitialize(self, problem):
        RheologyIncompressibleElasticity.preinitialize(self, problem)

        ModuleLinearElasticity.useReferenceState(self, self.useReferenceState)

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

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
