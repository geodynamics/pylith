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


class Solution(PetscComponent):
    """
    Abstract base class for solution field for problem.
    """

    import pythia.pyre.inventory

    from .SolnDisp import SolnDisp
    from .SolutionSubfield import subfieldFactory
    subfields = pythia.pyre.inventory.facilityArray("subfields", family="soln_subfields",
                                                    itemFactory=subfieldFactory, factory=SolnDisp)
    subfields.meta['tip'] = "Subfields in solution."

    def __init__(self, name="solution"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="solution")
        self.field = None

    def preinitialize(self, problem, mesh):
        """Do minimal initialization of solution.
        """
        from pylith.mpi.Communicator import mpi_is_root
        isRoot = mpi_is_root()
        if isRoot:
            self._info.log("Performing minimal initialization of solution.")

        from pylith.topology.Field import Field
        self.field = Field(mesh)
        self.field.setLabel("solution")
        spaceDim = mesh.getCoordSys().getSpaceDim()
        for subfield in self.subfields.components():
            subfield.initialize(problem.normalizer, spaceDim)
            if isRoot:
                self._debug.log("Adding subfield '%s' as '%s' with components %s to solution." %
                                (subfield.fieldName, subfield.userAlias, subfield.componentNames))
            descriptor = subfield.getTraitDescriptor("quadrature_order")
            if hasattr(descriptor.locator, "source") and descriptor.locator.source == "default":
                quadOrder = problem.defaults.quadOrder
            else:
                quadOrder = subfield.quadOrder
            self.field.subfieldAdd(subfield.fieldName, subfield.userAlias, subfield.vectorFieldType, 
                                   subfield.componentNames, subfield.scale.value, subfield.basisOrder, 
                                   quadOrder, subfield.dimension, subfield.isFaultOnly,
                                   subfield.cellBasis, subfield.feSpace, subfield.isBasisContinuous)

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)

    def _cleanup(self):
        if self.field:
            self.field.deallocate()

# FACTORIES ////////////////////////////////////////////////////////////


def solution():
    """Factory associated with Solution.
    """
    return Solution()


# End of file
