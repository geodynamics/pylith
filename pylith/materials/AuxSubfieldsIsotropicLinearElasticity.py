# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/materials/AuxSubfieldsIsotropicLinearElasticity.py
#
# @brief Python container for isotropic, linear elasticity subfields.

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsIsotropicLinearElasticity(PetscComponent):
    """
    Python container for isotropic, linear elasticity subfields.

    INVENTORY

    Properties
      - None

    Facilities
      - *shear_modulus* Shear modulus subfield.
      - *bulk_modulus* Bulk modulus subfield.
      - *reference_stress* Reference stress subfield.
      - *references_strain* Reference strain.
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    shearModulus = pythia.pyre.inventory.facility("shear_modulus", family="auxiliary_subfield", factory=Subfield)
    shearModulus.meta['tip'] = "Shear modulus subfield."

    bulkModulus = pythia.pyre.inventory.facility("bulk_modulus", family="auxiliary_subfield", factory=Subfield)
    bulkModulus.meta['tip'] = "Bulk modulus subfield."

    referenceStress = pythia.pyre.inventory.facility("reference_stress", family="auxiliary_subfield", factory=Subfield)
    referenceStress.meta['tip'] = "Reference stress subfield."

    referenceStrain = pythia.pyre.inventory.facility("reference_strain", family="auxiliary_subfield", factory=Subfield)
    referenceStrain.meta['tip'] = "Reference strain subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxsubfieldsisotropiclinearelasticity"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_subfields")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        PetscComponent._configure(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfields():
    """
    Factory associated with AuxSubfieldsIsotropicLinearElasticity.
    """
    return AuxSubfieldsIsotropicLinearElasticity()


# End of file
