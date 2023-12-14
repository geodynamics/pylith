# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.materials.DerivedSubfieldsElasticity import DerivedSubfieldsElasticity

class DerivedSubfieldsPoroelasticity(DerivedSubfieldsElasticity):

    """
    Derived subfields associated with the poroelasticity equation.

    For poroelastic materials these derived subfields are available for output in addition to the solution and auxiliary subfields.
    """
    DOC_CONFIG = {
        "cfg": """
            # The basis order for stress and strain should be at least 1 less than the basis order for displacement.
            [pylithapp.problem.materials.mat_elastic.derived_subfields]
            cauchy_stress.basis_order = 0
            cauchy_strain.basis_order = 0
            bulk_density.basis_order = 1
            water_content.basis_order = 1
        """
    }

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    bulkDensity = pythia.pyre.inventory.facility("bulk_density", family="subfield", factory=Subfield)
    bulkDensity.meta['tip'] = "bulk density subfield."

    waterContent = pythia.pyre.inventory.facility("water_content", family="subfield", factory=Subfield)
    waterContent.meta['tip'] = "water content subfield"


    def __init__(self, name="derivedsubfieldsporoelasticity"):
        """Constructor.
        """
        DerivedSubfieldsElasticity.__init__(self, name)

    def _defaults(self):
        self.bulkDensity.basis_order = 1
        self.waterContent.basis_order = 1

# FACTORIES ////////////////////////////////////////////////////////////

def derived_subfields():
    """Factory associated with DerivedSubfieldsPoroelasticity.
    """
    return DerivedSubfieldsPoroelasticity()


# End of file