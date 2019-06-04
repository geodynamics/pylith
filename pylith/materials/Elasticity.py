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
# @file pylith/materials/Elasticity.py
#
# @brief Python object for solving the elasticity equation.
#
# Factory: material

from .Material import Material
from .materials import Elasticity as ModuleElasticity

from .IsotropicLinearElasticity import IsotropicLinearElasticity


class Elasticity(Material, ModuleElasticity):
    """
    Python material property manager.

    INVENTORY

    Properties
      - *use_inertia* Include inertial term in elasticity equation.
      - *use_body_force* Include body force term in elasticity equation.

    Facilities
      - *bulk_rheology* Bulk rheology for elastic material.

    FACTORY: material
    """

    import pyre.inventory

    from pylith.topology.Subfield import subfieldFactory
    from pylith.utils.EmptyBin import EmptyBin

    from .AuxSubfieldsElasticity import AuxSubfieldsElasticity
    auxiliarySubfields = pyre.inventory.facilityArray(
        "auxiliary_subfields", itemFactory=subfieldFactory, factory=AuxSubfieldsElasticity)
    auxiliarySubfields.meta['tip'] = "Discretization of elasticity properties."

    from .DerivedSubfieldsElasticity import DerivedSubfieldsElasticity
    derivedSubfields = pyre.inventory.facilityArray(
        "derived_subfields", itemFactory=subfieldFactory, factory=DerivedSubfieldsElasticity)
    derivedSubfields.meta['tip'] = "Discretization of derived subfields (e.g., stress and strain)."

    useInertia = pyre.inventory.bool("use_inertia", default=False)
    useInertia.meta['tip'] = "Include inertial term in elasticity equation."

    useBodyForce = pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force term in elasticity equation."

    rheology = pyre.inventory.facility("bulk_rheology", family="elasticity_rheology", factory=IsotropicLinearElasticity)
    rheology.meta['tip'] = "Bulk rheology for elastic material."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="elasticity"):
        """
        Constructor.
        """
        Material.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Setup material.
        """
        self.rheology.preinitialize(mesh)
        Material.preinitialize(self, mesh)

        self.rheology.addAuxiliarySubfields(self)

        ModuleElasticity.useInertia(self, self.useInertia)
        ModuleElasticity.useBodyForce(self, self.useBodyForce)
        return

    def _createModuleObj(self):
        """
        Create handle to C++ Elasticity.
        """
        ModuleElasticity.__init__(self)
        ModuleElasticity.setBulkRheology(self, self.rheology)  # Material sets auxiliary db in rheology.
        return


# Factories

def material():
    """
    Factory associated with Elasticity.
    """
    return Elasticity()


# End of file
