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
# @file pylith/materials/IncompressibleElasticity.py
#
# @brief Python object for solving the elasticity equation with an incompressible material.
#
# Factory: material

from .Material import Material
from .materials import IncompressibleElasticity as ModuleIncompressibleElasticity


# VALIDATORS ///////////////////////////////////////////////////////////

class IncompressibleElasticity(Material, ModuleIncompressibleElasticity):
    """
    Python material property manager.

    INVENTORY

    Properties
      - *use_inertia* Include inertial term in elasticity equation.
      - *use_body_force* Include body force term in elasticity equation.

    Facilities
      - *bulk_rheology* Bulk rheology for incompessible elastic material.

    FACTORY: material
    """

    import pyre.inventory

    useInertia = pyre.inventory.bool("use_inertia", default=False)
    useInertia.meta['tip'] = "Include inertial term in elasticity equation."

    useBodyForce = pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force term in elasticity equation."

    rheology = pyre.inventory.facility(
        "bulk_rheology", familty="incompressible_elasticity_rheology", factory=IsotropicLinearElasticity)
    rheology.meta['tip'] = "Bulk rheology for elastic material."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="incompressibleelasticity"):
        """
        Constructor.
        """
        Material.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Setup material.
        """
        Material.preinitialize(self, mesh)

        ModuleIncompressibleElasticity.useIntertia(self, self.userInertia)
        ModuleIncompressibleElasticity.useBodyForce(self, self.useBodyForce)
        ModuleIncompressibleElasticity.setBulkRheology(self, self.rheology)

        return

# End of file
