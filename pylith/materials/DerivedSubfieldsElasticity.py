# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

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
