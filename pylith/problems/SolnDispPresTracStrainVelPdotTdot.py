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

# @file pylith/problems/SolnDispPresTracStrainVelTdotPdot.py
#
# @brief Python subfields container with displacement, pore pressure, and trace strain subfields, along with their time derivatives.

from pylith.utils.PetscComponent import PetscComponent
from .Solution import Solution as SolutionBase


class SolnDispPresTracStrainVelTdotPdot(PetscComponent):
    """Python subfields container with displacement, pore pressure, and trace strain subfields.

    IMPORTANT: Use the Solution class (below) to set this object as the default facilities array for the solution
    subfields.
    """

    import pythia.pyre.inventory

    from .SubfieldDisplacement import SubfieldDisplacement
    displacement = pythia.pyre.inventory.facility(
        "displacement", family="soln_subfield", factory=SubfieldDisplacement)
    displacement.meta['tip'] = "Displacement subfield."

    from .SubfieldPressure import SubfieldPressure
    pressure = pythia.pyre.inventory.facility(
        "pressure", family="soln_subfield", factory=SubfieldPressure)
    pressure.meta['tip'] = "Pressure subfield."

    from .SubfieldTraceStrain import SubfieldTraceStrain
    trace_strain = pythia.pyre.inventory.facility(
        "trace_strain", family="soln_subfield", factory=SubfieldTraceStrain)
    trace_strain.meta['tip'] = "Trace strain subfield."

    from .SubfieldVelocity import SubfieldVelocity
    velocity = pythia.pyre.inventory.facility(
        "velocity", family="soln_subfield", factory=SubfieldVelocity)
    velocity.meta['tip'] = "Velocity subfield."

    from .SubfieldPressureDot import SubfieldPressureDot
    pressure_t = pythia.pyre.inventory.facility(
        "pressure_t", family="soln_subfield", factory=SubfieldPressureDot)
    pressure_t.meta['tip'] = "Pressure_t subfield."

    from .SubfieldTraceStrainDot import SubfieldTraceStrainDot
    trace_strain_t = pythia.pyre.inventory.facility(
        "trace_strain_t", family="soln_subfield", factory=SubfieldTraceStrainDot)
    trace_strain_t.meta['tip'] = "Trace strain_t subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="solndispprestracstrainveltdotpdot"):
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
        return [self.displacement, self.pressure, self.trace_strain, self.velocity, self.pressure_t, self.trace_strain_t]


class Solution(SolutionBase):
    """Python solution field with displacement, pressure, and trace strain subfields.
    """

    import pythia.pyre.inventory

    from .SolutionSubfield import subfieldFactory
    subfields = pythia.pyre.inventory.facilityArray(
        "subfields", family="soln_subfields", itemFactory=subfieldFactory, factory=SolnDispPresTracStrainVelTdotPdot)
    subfields.meta['tip'] = "Subfields in solution."


# FACTORIES ////////////////////////////////////////////////////////////
def solution():
    """Factory associated with Solution.
    """
    return Solution()


# End of file
