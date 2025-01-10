# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .Material import Material
from .materials import Elasticity as ModuleElasticity

from .IsotropicLinearElasticity import IsotropicLinearElasticity


class Elasticity(Material, ModuleElasticity):
    """
    Material behavior governed by the elasticity equation.

    Implements `Material`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_elastic]
            label_value = 4
            description = Upper crust elastic material
            use_body_force = False
            bulk_rheology = pylith.materials.IsotropicLinearElasticity

            auxiliary_subfields.density.basis_order = 0
            derived_subfields.cauchy_stress.basis_order = 1
            derived_subfields.cauchy_strain.basis_order = 1
        """
    }

    import pythia.pyre.inventory

    useBodyForce = pythia.pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force term in elasticity equation."

    rheology = pythia.pyre.inventory.facility("bulk_rheology", family="elasticity_rheology", factory=IsotropicLinearElasticity)
    rheology.meta['tip'] = "Bulk rheology for elastic material."

    def __init__(self, name="elasticity"):
        """Constructor.
        """
        Material.__init__(self, name)

    def _defaults(self):
        from .AuxSubfieldsElasticity import AuxSubfieldsElasticity
        self.auxiliarySubfields = AuxSubfieldsElasticity("auxiliary_subfields")

        from .DerivedSubfieldsElasticity import DerivedSubfieldsElasticity
        self.derivedSubfields = DerivedSubfieldsElasticity("derived_subfields")

    def preinitialize(self, problem):
        """Setup material.
        """
        self.rheology.preinitialize(problem)
        Material.preinitialize(self, problem)
        self.rheology.addAuxiliarySubfields(self, problem)
        ModuleElasticity.useBodyForce(self, self.useBodyForce)

    def _createModuleObj(self):
        """Create handle to C++ Elasticity.
        """
        ModuleElasticity.__init__(self)
        ModuleElasticity.setBulkRheology(self, self.rheology)  # Material sets auxiliary db in rheology.


# Factories

def material():
    """Factory associated with Elasticity.
    """
    return Elasticity()


# End of file
