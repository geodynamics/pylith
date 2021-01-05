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
# @file pylith/materials/RheologyElasticity.py
#
# @brief Python material for isotropic, linearly elastic, plane
# strain material.
#
# Factory: elasticity_rheology

from pylith.utils.PetscComponent import PetscComponent
from .materials import RheologyElasticity as ModuleRheology


class RheologyElasticity(PetscComponent, ModuleRheology):
    """
    Python object for bulk rheology of an elastic material.

    INVENTORY

    Properties
      - None

    Facilities
      - *auxiliary_subfields* Discretization of physical properties and state variables.

    FACTORY: elasticity_rheology
    """
    import pyre.inventory

    from pylith.topology.Subfield import subfieldFactory
    from pylith.utils.EmptyBin import EmptyBin

    auxiliarySubfields = pyre.inventory.facilityArray(
        "auxiliary_subfields", itemFactory=subfieldFactory, factory=EmptyBin)
    auxiliarySubfields.meta['tip'] = "Discretization information for physical properties and state variables."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="rheologyelasticity")
        return

    def preinitialize(self, problem):
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Performing minimal initialization of elasticity rheology '%s'." % self.aliases[-1])

        self._createModuleObj()
        return

    def addAuxiliarySubfields(self, material, problem):
        for subfield in self.auxiliarySubfields.components():
            fieldName = subfield.aliases[-1]
            quadOrder = problem.defaults.quadOrder if subfield.quadOrder < 0 else subfield.quadOrder
            material.setAuxiliarySubfieldDiscretization(fieldName, subfield.basisOrder, quadOrder, subfield.dimension,
                                                        subfield.cellBasis, subfield.isBasisContinuous, subfield.feSpace)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        raise NotImplementedError("Implement in derived class.")


# FACTORIES ////////////////////////////////////////////////////////////

def elasticity_rheology():
    """
    Factory associated with RheologyElasticity.
    """
    return RheologyElasticity()


# End of file
