# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

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
        from pylith.mpi.Communicator import mpi_is_root
        if mpi_is_root():
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
