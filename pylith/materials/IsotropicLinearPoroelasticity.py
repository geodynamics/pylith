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
# @file pylith/materials/IsotropicLinearPoroelasticity.py
#
# @brief Python material for isotropic, linearly elastic, plane
# strain material.
#
# Factory: poroelasticity_rheology

from .RheologyPoroelasticity import RheologyPoroelasticity
from .materials import IsotropicLinearPoroelasticity as ModuleLinearPoroelasticity


class IsotropicLinearPoroelasticity(RheologyPoroelasticity, ModuleLinearPoroelasticity):
    """Python material for isotropic, linearly poroelastic plane strain.

    FACTORY: poroelasticity_rheology
    """

    import pythia.pyre.inventory

    useReferenceState = pythia.pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    useTensorPermeability = pythia.pyre.inventory.bool("use_tensor_permeability", default=False)
    useTensorPermeability.meta['tip'] = "Use tensor permeability."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="isotropiclinearporoelasticity"):
        """Constructor.
        """
        RheologyPoroelasticity.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsIsotropicLinearPoroelasticity import AuxSubfieldsIsotropicLinearPoroelasticity
        self.auxiliarySubfields = AuxSubfieldsIsotropicLinearPoroelasticity("auxiliary_subfields")

    def preinitialize(self, mesh):
        RheologyPoroelasticity.preinitialize(self, mesh)

        ModuleLinearPoroelasticity.useReferenceState(self, self.useReferenceState)
        ModuleLinearPoroelasticity.useTensorPermeability(self, self.useTensorPermeability)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleLinearPoroelasticity.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def poroelasticity_rheology():
    """Factory associated with IsotropicLinearPoroelasticity.
    """
    return IsotropicLinearPoroelasticity()


# End of file
