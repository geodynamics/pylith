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


class SubfieldTraceStrain(SolutionSubfield):
    """
    Object for defining attributes of the trace strain solution subfield.

    Implements `SolutionSubfield`.
    """
    DOC_CONFIG = {
        "cfg": """
        [pylithapp.problems.solution.subfields.trace_strain]
        alias = trace_strain
        basis_order = 1
        """
    }

    fieldName = "trace_strain"

    def __init__(self, name="subfieldtracestrain"):
        """Constructor.
        """
        SolutionSubfield.__init__(self, name)
        return

    def _defaults(self):
        self.userAlias = self.fieldName

    def initialize(self, normalizer, spaceDim):
        """Initialize subfield metadata.
        """
        from pylith.topology.Field import Field
        from pythia.pyre.units.unit import one
        self.vectorFieldType = Field.SCALAR
        self.scale = one
        self._setComponents(spaceDim)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        SolutionSubfield._configure(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def soln_subfield():
    """Factory associated with SubfieldTraceStrain.
    """
    return SubfieldTraceStrain()


# End of file
