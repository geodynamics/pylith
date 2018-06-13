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
# @file pylith/materials/IsotropicLinearMaxwellPlaneStrain.py
#
# @brief Python material for isotropic, linearly Maxwell viscoelastic, plane
# strain material.
#
# Factory: material

from .Material import Material
from .materials import IsotropicLinearMaxwellPlaneStrain as ModuleMaterial


class IsotropicLinearMaxwellPlaneStrain(Material, ModuleMaterial):
    """
    Python material for isotropic, linear Maxwell viscoelastic plane strain.

    INVENTORY

    Properties
      - *use_inertia* Include inertial term in elasticity equation.
      - *use_body_force* Include body force term in elasticity equation.
      - *use_reference_state* Use reference stress/strain state.

    Facilities
      - *auxiliary_subfields* Discretization of physical properties and state variables.

    FACTORY: material
    """

    import pyre.inventory

    useInertia = pyre.inventory.bool("use_inertia", default=False)
    useInertia.meta['tip'] = "Include inertial term in elasticity equation."

    useBodyForce = pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force term in elasticity equation."

    useReferenceState = pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    from .AuxFieldsIsotropicLinearMaxwell import AuxFieldsIsotropicLinearMaxwell
    from pylith.topology.AuxSubfield import subfieldFactory
    auxSubfields = pyre.inventory.facilityArray("auxiliary_subfields", itemFactory=subfieldFactory, factory=AuxFieldsIsotropicLinearMaxwell)
    auxSubfields.meta['tip'] = "Discretization of physical properties and state variables."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="isotropiclinearmaxwellplanestrain"):
        """
        Constructor.
        """
        Material.__init__(self, name)
        return

    def preinitialize(self, mesh):
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Performing minimal initialization of material '%s'" % self.aliases[-1])

        Material.preinitialize(self, mesh)

        ModuleMaterial.useInertia(self, self.useInertia)
        ModuleMaterial.useBodyForce(self, self.useBodyForce)
        ModuleMaterial.useReferenceState(self, self.useReferenceState)

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        Material._configure(self)
        return

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        ModuleMaterial.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def material():
    """
    Factory associated with IsotropicLinearMaxwellPlaneStrain.
    """
    return IsotropicLinearMaxwellPlaneStrain()


# End of file
