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
# @file pylith/materials/Poroelasticity.py
#
# @brief Python object for solving the Poroelasticity equation.
#
# Factory: material

from .Material import Material
from .materials import Poroelasticity as ModulePoroelasticity

from .IsotropicLinearPoroelasticity import IsotropicLinearPoroelasticity


class Poroelasticity(Material, ModulePoroelasticity):
    """
    Python material property manager.

    INVENTORY

    Properties
      - *use_inertia* Include inertial term in Poroelasticity equation.
      - *use_body_force* Include body force term in Poroelasticity equation.

    Facilities
      - *bulk_rheology* Bulk rheology for poroelastic material.

    FACTORY: material
    """

    import pythia.pyre.inventory

    useInertia = pythia.pyre.inventory.bool("use_inertia", default=False)
    useInertia.meta['tip'] = "Include inertia term in Poroelasticity equation."

    useBodyForce = pythia.pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force term in Poroelasticity equation."

    useSourceDensity = pythia.pyre.inventory.bool("use_source_density", default=False)
    useSourceDensity.meta['tip'] = "Include source density term in Poroelasticity equation."

    rheology = pythia.pyre.inventory.facility(
        "bulk_rheology", family="poroelasticity_rheology", factory=IsotropicLinearPoroelasticity)
    rheology.meta['tip'] = "Bulk rheology for poroelastic material."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="poroelasticity"):
        """Constructor.
        """
        Material.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsPoroelasticity import AuxSubfieldsPoroelasticity
        self.auxiliarySubfields = AuxSubfieldsPoroelasticity("auxiliary_subfields")

        from .DerivedSubfieldsElasticity import DerivedSubfieldsElasticity
        self.derivedSubfields = DerivedSubfieldsElasticity("derived_subfields")

    def preinitialize(self, problem):
        """Setup material.
        """
        self.rheology.preinitialize(problem)
        Material.preinitialize(self, problem)

        self.rheology.addAuxiliarySubfields(self, problem)

        ModulePoroelasticity.useInertia(self, self.useInertia)
        ModulePoroelasticity.useBodyForce(self, self.useBodyForce)
        ModulePoroelasticity.useSourceDensity(self, self.useSourceDensity)
        return

    def _createModuleObj(self):
        """Create handle to C++ Poroelasticity.
        """
        ModulePoroelasticity.__init__(self)
        ModulePoroelasticity.setBulkRheology(self, self.rheology)  # Material sets auxiliary db in rheology.
        return


# Factories

def material():
    """Factory associated with Poroelasticity.
    """
    return Poroelasticity()


# End of file
