# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsIsotropicPowerLaw(PetscComponent):
    """
    Auxiliary subfields associated with the isotropic power-law viscoelastic bulk rheology.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_powerlaw.rheology.auxiliary_fields]
            shear_modulus.basis_order = 1
            bulk_modulus.basis_order = 1
            power_law_reference_strain_rate = 1
            power_law_reference_stress = 1
            power_law_exponent.basis_order = 1
            viscous_strain.basis_order = 1
            deviatoric_stress.basis_order = 1
            reference_stress.basis_order = 0
            reference_strain.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    shearModulus = pythia.pyre.inventory.facility("shear_modulus", family="auxiliary_subfield", factory=Subfield)
    shearModulus.meta['tip'] = "Shear modulus subfield."

    bulkModulus = pythia.pyre.inventory.facility("bulk_modulus", family="auxiliary_subfield", factory=Subfield)
    bulkModulus.meta['tip'] = "Bulk modulus subfield."

    powerLawReferenceStrainRate = pythia.pyre.inventory.facility("power_law_reference_strain_rate", family="auxiliary_subfield", factory=Subfield)
    powerLawReferenceStrainRate.meta['tip'] = "Power-law reference strain rate subfield."

    powerLawReferenceStress = pythia.pyre.inventory.facility("power_law_reference_stress", family="auxiliary_subfield", factory=Subfield)
    powerLawReferenceStress.meta['tip'] = "Power-law reference stress subfield."

    powerLawExponent = pythia.pyre.inventory.facility("power_law_exponent", family="auxiliary_subfield", factory=Subfield)
    powerLawExponent.meta['tip'] = "Power-law exponent subfield."

    viscousStrain = pythia.pyre.inventory.facility("viscous_strain", family="auxiliary_subfield", factory=Subfield)
    viscousStrain.meta['tip'] = "Viscous strain subfield."

    deviatoricStress = pythia.pyre.inventory.facility("deviatoric_stress", family="auxiliary_subfield", factory=Subfield)
    deviatoricStress.meta['tip'] = "Deviatoric stress subfield."

    referenceStress = pythia.pyre.inventory.facility("reference_stress", family="auxiliary_subfield", factory=Subfield)
    referenceStress.meta['tip'] = "Reference stress subfield."

    referenceStrain = pythia.pyre.inventory.facility("reference_strain", family="auxiliary_subfield", factory=Subfield)
    referenceStrain.meta['tip'] = "Reference strain subfield."

    def __init__(self, name="auxfieldsisotropicpowerlaw"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_subfields")

    def _configure(self):
        PetscComponent._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfields():
    """Factory associated with AuxSubfieldsIsotropicPowerLaw.
    """
    return AuxSubfieldsIsotropicPowerLaw()


# End of file
