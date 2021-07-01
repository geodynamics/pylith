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

# @file pylith/materials/AuxSubfieldsIsotropicPowerLaw.py
##
# @brief Python subfields container for isotropic power-law
# viscoelastic subfields.

from pylith.utils.PetscComponent import PetscComponent

# AuxSubfieldsIsotropicPowerLaw class


class AuxSubfieldsIsotropicPowerLaw(PetscComponent):
    """Python container for isotropic power-law viscoelastic subfields.

    FACTORY: auxiliary_subfields
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    shearModulus = pythia.pyre.inventory.facility("shear_modulus", family="auxiliary_subfield", factory=Subfield)
    shearModulus.meta['tip'] = "Shear modulus subfield."

    bulkModulus = pythia.pyre.inventory.facility("bulk_modulus", family="auxiliary_subfield", factory=Subfield)
    bulkModulus.meta['tip'] = "Bulk modulus subfield."

    powerLawReferenceStrainRate = pythia.pyre.inventory.facility("power_law_reference_strain_rate", family="auxiliary_subfield",
                                                                 factory=Subfield)
    powerLawReferenceStrainRate.meta['tip'] = "Power-law reference strain rate subfield."

    powerLawReferenceStress = pythia.pyre.inventory.facility("power_law_reference_stress", family="auxiliary_subfield",
                                                             factory=Subfield)
    powerLawReferenceStress.meta['tip'] = "Power-law reference stress subfield."

    powerLawExponent = pythia.pyre.inventory.facility("power_law_exponent", family="auxiliary_subfield", factory=Subfield)
    powerLawExponent.meta['tip'] = "Power-law exponent subfield."

    viscousStrain = pythia.pyre.inventory.facility("viscous_strain", family="auxiliary_subfield", factory=Subfield)
    viscousStrain.meta['tip'] = "Viscous strain subfield."

    stress = pythia.pyre.inventory.facility("stress", family="auxiliary_subfield", factory=Subfield)
    stress.meta['tip'] = "Stress subfield."

    referenceStress = pythia.pyre.inventory.facility("reference_stress", family="auxiliary_subfield", factory=Subfield)
    referenceStress.meta['tip'] = "Reference stress subfield."

    referenceStrain = pythia.pyre.inventory.facility("reference_strain", family="auxiliary_subfield", factory=Subfield)
    referenceStrain.meta['tip'] = "Reference strain subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxfieldsisotropicpowerlaw"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_subfields")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        PetscComponent._configure(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfields():
    """Factory associated with AuxSubfieldsIsotropicPowerLaw.
    """
    return AuxSubfieldsIsotropicPowerLaw()


# End of file
