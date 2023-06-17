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
#

# @file pylith/problems/SolnDispPresTracStrainVelPdotTdotLagrange.py
#
# @brief Python subfields container with displacement, pore pressure, and trace strain subfields, along with their time derivatives.

from pylith.utils.PetscComponent import PetscComponent
from .Solution import Solution as SolutionBase


class SolnDispPresTracStrainVelPdotTdotLagrange(PetscComponent):
    """
        Python solution field with displacement, pressure, and trace strain subfields, along with their time derivatives, and a lagrange fault.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem]
            solution = pylith.problems.SolnDispPresTracStrainVelPdotTdotLagrange
        """
    }

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
    traceStrain = pythia.pyre.inventory.facility(
        "trace_strain", family="soln_subfield", factory=SubfieldTraceStrain)
    traceStrain.meta['tip'] = "Trace strain subfield."

    from .SubfieldVelocity import SubfieldVelocity
    velocity = pythia.pyre.inventory.facility(
        "velocity", family="soln_subfield", factory=SubfieldVelocity)
    velocity.meta['tip'] = "Velocity subfield."

    from .SubfieldPressureDot import SubfieldPressureDot
    pressureT = pythia.pyre.inventory.facility(
        "pressure_t", family="soln_subfield", factory=SubfieldPressureDot)
    pressureT.meta['tip'] = "PressureT subfield."

    from .SubfieldTraceStrainDot import SubfieldTraceStrainDot
    traceStrainT = pythia.pyre.inventory.facility(
        "trace_strain_t", family="soln_subfield", factory=SubfieldTraceStrainDot)
    traceStrainT.meta['tip'] = "TraceStrainT subfield."

    from .SubfieldLagrangeFault import SubfieldLagrangeFault
    lagrangeFault = pythia.pyre.inventory.facility("lagrange_fault", family="soln_subfield", factory=SubfieldLagrangeFault)
    lagrangeFault.meta['tip'] = "Fault Lagrange multiplier subfield."

    def __init__(self, name="solndispprestracstrainvelpdottdotlagrange"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="soln_subfields")

    def _configure(self):
        PetscComponent._configure(self)

    def components(self):
        """Order of facilities in Inventory is ambiguous, so overwrite
        components() to insure order is [displacement, pressure, trace_strain, velocity, pressure_t, trace_strain_t, lagrange_fault].

        """
        return [self.displacement, self.pressure, self.trace_strain, self.velocity, self.pressureT, self.traceStrainT, self.lagrangeFault]


class Solution(SolutionBase):
    """Python solution field with displacement, pressure, and trace strain subfields, along with their time derivatives, and a lagrange fault.
    """

    import pythia.pyre.inventory

    from .SolutionSubfield import subfieldFactory
    subfields = pythia.pyre.inventory.facilityArray(
        "subfields", family="soln_subfields", itemFactory=subfieldFactory, factory=SolnDispPresTracStrainVelPdotTdotLagrange)
    subfields.meta['tip'] = "Subfields in solution."


# FACTORIES ////////////////////////////////////////////////////////////
def solution():
    """Factory associated with Solution.
    """
    return Solution()


# End of file