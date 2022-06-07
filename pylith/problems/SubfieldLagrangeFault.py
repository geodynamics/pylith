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

from .SolutionSubfield import SolutionSubfield


class SubfieldLagrangeFault(SolutionSubfield):
    """
    Object for defining attributes of the fault Lagrange multiplier solution subfield.

    Implements `SolutionSubfield`.
    """
    DOC_CONFIG = {
        "cfg": """
        [pylithapp.problems.solution.subfields.lagrange_fault]
        alias = lagrange_multiplier_fault
        basis_order = 1
        """
    }

    fieldName = "lagrange_multiplier_fault"

    def __init__(self, name="subfieldlagrangefault"):
        """Constructor.
        """
        SolutionSubfield.__init__(self, name)

    def _defaults(self):
        self.userAlias = self.fieldName

    def initialize(self, normalizer, spaceDim):
        """Initialize subfield metadata.
        """
        from pylith.topology.Field import Field
        self.dimension = spaceDim - 1
        self.vectorFieldType = Field.VECTOR
        self.scale = normalizer.getPressureScale()
        self._setComponents(spaceDim)
        self.isFaultOnly = True

    def _configure(self):
        """Set members based using inventory.
        """
        SolutionSubfield._configure(self)

# FACTORIES ////////////////////////////////////////////////////////////


def soln_subfield():
    """Factory associated with SubfieldLagrangeFault.
    """
    return SubfieldLagrangeFault()


# End of file
