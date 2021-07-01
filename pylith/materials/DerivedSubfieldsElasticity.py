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
# @file pylith/materials/AuxSubieldsElasticity.py
#
# @brief Python container for elasticity equation subfields.

from pylith.utils.PetscComponent import PetscComponent


class DerivedSubfieldsElasticity(PetscComponent):
    """Python container for derived subfields for elasticity.

    FACTORY: derived_subfields
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    cauchyStress = pythia.pyre.inventory.facility("cauchy_stress", family="subfield", factory=Subfield)
    cauchyStress.meta['tip'] = "Cauchy stress subfield."

    cauchyStrain = pythia.pyre.inventory.facility("cauchy_strain", family="subfield", factory=Subfield)
    cauchyStrain.meta['tip'] = "Cauchy strain subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="derivedsubfieldselasticity"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="derived_subfields")
        return


# FACTORIES ////////////////////////////////////////////////////////////

def derived_subfields():
    """Factory associated with DerivedSubfieldsElasticity.
    """
    return DerivedSubfieldsElasticity()


# End of file
