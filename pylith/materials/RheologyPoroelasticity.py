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

from pylith.utils.PetscComponent import PetscComponent
from .materials import RheologyPoroelasticity as ModuleRheology


class RheologyPoroelasticity(PetscComponent, ModuleRheology):
    """
    Abstract base class for bulk rheology of poroelastic material.
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import subfieldFactory
    from pylith.utils.EmptyBin import EmptyBin

    auxiliarySubfields = pythia.pyre.inventory.facilityArray("auxiliary_subfields", itemFactory=subfieldFactory, factory=EmptyBin)
    auxiliarySubfields.meta['tip'] = "Discretization information for physical properties and state variables."

    def __init__(self, name="rheologyporoelasticity"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="rheologyporoelasticity")

    def preinitialize(self, problem):
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Performing minimal initialization of poroelasticity rheology '%s'." %
                           self.aliases[-1])

        self._createModuleObj()

    def addAuxiliarySubfields(self, material, problem):
        for subfield in self.auxiliarySubfields.components():
            fieldName = subfield.aliases[-1]
            descriptor = subfield.getTraitDescriptor("quadrature_order")
            if hasattr(descriptor.locator, "source") and descriptor.locator.source == "default":
                quadOrder = problem.defaults.quadOrder
            else:
                quadOrder = subfield.quadOrder
            material.setAuxiliarySubfieldDiscretization(fieldName, subfield.basisOrder, quadOrder, subfield.dimension,
                                                        subfield.cellBasis, subfield.feSpace, subfield.isBasisContinuous)

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        raise NotImplementedError("Implement in derived class.")


# End of file
