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


class DerivedSubfieldsElasticity(PetscComponent):
    """
    Derived subfields associated with the elasticity equation.

    For elastic materials these derived subfields are available for output in addition to the solution and auxiliary subfields.
    """
    DOC_CONFIG = {
        "cfg": """
            # The basis order for stress and strain should be at least 1 less than the basis order for displacement.
            [pylithapp.problem.materials.mat_elastic.derived_subfields]
            cauchy_stress.basis_order = 0
            cauchy_strain.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    cauchyStress = pythia.pyre.inventory.facility("cauchy_stress", family="subfield", factory=Subfield)
    cauchyStress.meta['tip'] = "Cauchy stress subfield."

    cauchyStrain = pythia.pyre.inventory.facility("cauchy_strain", family="subfield", factory=Subfield)
    cauchyStrain.meta['tip'] = "Cauchy strain subfield."

    def __init__(self, name="derivedsubfieldselasticity"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="derived_subfields")

    def _defaults(self):
        self.cauchyStress.basisOrder = 0
        self.cauchyStrain.basisOrder = 0


# FACTORIES ////////////////////////////////////////////////////////////

def derived_subfields():
    """Factory associated with DerivedSubfieldsElasticity.
    """
    return DerivedSubfieldsElasticity()


# End of file
