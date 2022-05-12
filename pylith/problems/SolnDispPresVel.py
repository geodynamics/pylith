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
from .Solution import Solution as SolutionBase


class SolnDispPresVel(PetscComponent):
    """
    Container for solution subfields with displacement, pore pressure, and velocity subfields.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem]
            solution = pylith.problems.SolnDispPresVel
        """
    }

    import pythia.pyre.inventory

    from .SubfieldDisplacement import SubfieldDisplacement
    displacement = pythia.pyre.inventory.facility("displacement", family="soln_subfield", factory=SubfieldDisplacement)
    displacement.meta['tip'] = "Displacement subfield."

    from .SubfieldPressure import SubfieldPressure
    pressure = pythia.pyre.inventory.facility("pressure", family="soln_subfield", factory=SubfieldPressure)
    pressure.meta['tip'] = "Pressure subfield."

    from .SubfieldVelocity import SubfieldVelocity
    velocity = pythia.pyre.inventory.facility("velocity", family="soln_subfield", factory=SubfieldVelocity)
    velocity.meta['tip'] = "Velocity subfield."

    def __init__(self, name="solndisppresvel"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="soln_subfields")

    def _configure(self):
        PetscComponent._configure(self)

    def components(self):
        """Order of facilities in Inventory is ambiguous, so overwrite
        components() to insure order is [displacement, pressure, velocity].
        """
        return [self.displacement, self.pressure, self.velocity]


class Solution(SolutionBase):
    """Python solution field with displacement, pressure, and trace strain subfields.
    """

    import pythia.pyre.inventory

    from .SolutionSubfield import subfieldFactory
    subfields = pythia.pyre.inventory.facilityArray("subfields", family="soln_subfields", itemFactory=subfieldFactory, factory=SolnDispPresVel)
    subfields.meta['tip'] = "Subfields in solution."


# FACTORIES ////////////////////////////////////////////////////////////
def solution():
    """Factory associated with Solution.
    """
    return Solution()


# End of file
