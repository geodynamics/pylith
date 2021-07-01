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

# @file pylith/problems/SolnDispPresTracStrain.py
#
# @brief Python subfields container with displacement, pore pressure, and trace strain subfields.

from pylith.utils.PetscComponent import PetscComponent
from .Solution import Solution as SolutionBase


class SolnDispPresTracStrain(PetscComponent):
    """Python subfields container with displacement, pore pressure, and trace strain subfields.

    IMPORTANT: Use the Solution class (below) to set this object as the default facilities array for the solution
    subfields.
    """

    import pythia.pyre.inventory

    from .SubfieldDisplacement import SubfieldDisplacement
    displacement = pythia.pyre.inventory.facility("displacement", family="soln_subfield", factory=SubfieldDisplacement)
    displacement.meta['tip'] = "Displacement subfield."

    from .SubfieldPressure import SubfieldPressure
    pressure = pythia.pyre.inventory.facility("pressure", family="soln_subfield", factory=SubfieldPressure)
    pressure.meta['tip'] = "Pressure subfield."

    from .SubfieldTraceStrain import SubfieldTraceStrain
    traceStrain = pythia.pyre.inventory.facility("trace_strain", family="soln_subfield", factory=SubfieldTraceStrain)
    traceStrain.meta['tip'] = "Trace strain subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="solndispprestracstrain"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="soln_subfields")
        return

    def _configure(self):
        PetscComponent._configure(self)
        return

    def components(self):
        """Order of facilities in Inventory is ambiguous, so overwrite
        components() to insure order is [displacement, pressure, trace_strain].

        """
        return [self.displacement, self.pressure, self.traceStrain]


class Solution(SolutionBase):
    """Python solution field with displacement, pressure, and trace strain subfields.
    """

    import pythia.pyre.inventory

    from .SolutionSubfield import subfieldFactory
    subfields = pythia.pyre.inventory.facilityArray(
        "subfields", family="soln_subfields", itemFactory=subfieldFactory, factory=SolnDispPresTracStrain)
    subfields.meta['tip'] = "Subfields in solution."


# FACTORIES ////////////////////////////////////////////////////////////
def solution():
    """Factory associated with Solution.
    """
    return Solution()


# End of file
