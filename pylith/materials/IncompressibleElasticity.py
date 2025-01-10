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
from .materials import IncompressibleElasticity as ModuleIncompressibleElasticity

from .IsotropicLinearIncompElasticity import IsotropicLinearIncompElasticity


class IncompressibleElasticity(Material, ModuleIncompressibleElasticity):
    """
    Material behavior governed by the elasticity equation.

    Implements `Material`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_incompelastic]
            description = Upper crust incompressible elastic material
            label_value = 3
            use_body_force = True
            bulk_rheology = pylith.materials.IsotropicLinearIncompElasticity

            auxiliary_subfields.density.basis_order = 0
            auxiliary_subfields.body_force.basis_order = 0
            derived_subfields.cauchy_stress.basis_order = 1
            derived_subfields.cauchy_strain.basis_order = 1
        """
    }

    import pythia.pyre.inventory

    useBodyForce = pythia.pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force term in elasticity equation."

    rheology = pythia.pyre.inventory.facility("bulk_rheology", family="incompressible_elasticity_rheology", factory=IsotropicLinearIncompElasticity)
    rheology.meta['tip'] = "Bulk rheology for elastic material."

    def __init__(self, name="incompressibleelasticity"):
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
        ModuleIncompressibleElasticity.useBodyForce(self, self.useBodyForce)

    def _createModuleObj(self):
        """Create handle to C++ IncompressibleElasticity.
        """
        ModuleIncompressibleElasticity.__init__(self)
        ModuleIncompressibleElasticity.setBulkRheology(self, self.rheology)  # Material sets auxiliary db in rheology.


# Factories

def material():
    """Factory associated with IncompressibleElasticity.
    """
    return IncompressibleElasticity()


# End of file
