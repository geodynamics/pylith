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

from .Material import Material
from .materials import Poroelasticity as ModulePoroelasticity

from .IsotropicLinearPoroelasticity import IsotropicLinearPoroelasticity


class Poroelasticity(Material, ModulePoroelasticity):
    """
    Material behavior governed by the poroelasticity equation.

    Implements `Material`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_poroelastic]
            description = Upper crust poroelastic material
            label_value = 3
            use_body_force = True
            use_source_density = False
            use_state_variables = True
            bulk_rheology = pylith.materials.IsotropicLinearPoroelasticity

            auxiliary_subfields.density.basis_order = 0
            auxiliary_subfields.body_force.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    useBodyForce = pythia.pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force term in Poroelasticity equation."

    useSourceDensity = pythia.pyre.inventory.bool("use_source_density", default=False)
    useSourceDensity.meta['tip'] = "Include source_density term in Poroelasticity equation."

    useStateVars = pythia.pyre.inventory.bool("use_state_variables", default=False)
    useStateVars.meta['tip'] = "Update porosity state variable using compaction formulation."

    rheology = pythia.pyre.inventory.facility("bulk_rheology", family="poroelasticity_rheology", factory=IsotropicLinearPoroelasticity)
    rheology.meta['tip'] = "Bulk rheology for poroelastic material."

    def __init__(self, name="poroelasticity"):
        """Constructor.
        """
        Material.__init__(self, name)

    def _defaults(self):
        from .AuxSubfieldsPoroelasticity import AuxSubfieldsPoroelasticity
        self.auxiliarySubfields = AuxSubfieldsPoroelasticity("auxiliary_subfields")

        from .DerivedSubfieldsPoroelasticity import DerivedSubfieldsPoroelasticity
        self.derivedSubfields = DerivedSubfieldsPoroelasticity("derived_subfields")

    def preinitialize(self, problem):
        """Setup material.
        """
        self.rheology.preinitialize(problem)
        Material.preinitialize(self, problem)

        self.rheology.addAuxiliarySubfields(self, problem)

        ModulePoroelasticity.useBodyForce(self, self.useBodyForce)
        ModulePoroelasticity.useSourceDensity(self, self.useSourceDensity)
        ModulePoroelasticity.useStateVars(self, self.useStateVars)        

    def _createModuleObj(self):
        """Create handle to C++ Poroelasticity.
        """
        ModulePoroelasticity.__init__(self)
        ModulePoroelasticity.setBulkRheology(self, self.rheology)  # Material sets auxiliary db in rheology.


# Factories

def material():
    """Factory associated with Poroelasticity.
    """
    return Poroelasticity()


# End of file
