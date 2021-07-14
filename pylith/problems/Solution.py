# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/problems/Solution.py
#
# @brief Python solution field for problem.
#
# Factory: solution.

from pylith.utils.PetscComponent import PetscComponent


class Solution(PetscComponent):
    """Python abstract base class for solution field for problem.

    FACTORY: solution.
    """

    import pythia.pyre.inventory

    from .SolnDisp import SolnDisp
    from .SolutionSubfield import subfieldFactory
    subfields = pythia.pyre.inventory.facilityArray("subfields", family="soln_subfields",
                                                    itemFactory=subfieldFactory, factory=SolnDisp)
    subfields.meta['tip'] = "Subfields in solution."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="solution"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="solution")
        self.field = None
        return

    def preinitialize(self, problem, mesh):
        """Do minimal initialization of solution.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Performing minimal initialization of solution.")

        from pylith.topology.Field import Field
        self.field = Field(mesh)
        self.field.setLabel("solution")
        spaceDim = mesh.getCoordSys().getSpaceDim()
        for subfield in self.subfields.components():
            subfield.initialize(problem.normalizer, spaceDim)
            if 0 == comm.rank:
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
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)
        return

    def _cleanup(self):
        if self.field:
            self.field.deallocate()
        return

# FACTORIES ////////////////////////////////////////////////////////////


def solution():
    """Factory associated with Solution.
    """
    return Solution()


# End of file
