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


class SubfieldTraceStrainDot(SolutionSubfield):
    """
    Object for defining attributes of the time derivative of trace strain solution subfield.

    Implements `SolutionSubfield`.
    """
    DOC_CONFIG = {
        "cfg": """
        [pylithapp.problems.solution.subfields.trace_strain_t]
        alias = trace_strain_t
        basis_order = 1
        """
    }

    fieldName = "trace_strain_t"

    def __init__(self, name="subfieldtracestrain_t"):
        """
        Constructor.
        """
        SolutionSubfield.__init__(self, name)

    def _defaults(self):
        self.userAlias = self.fieldName

    def initialize(self, normalizer, spaceDim):
        """
        Initialize subfield metadata.
        """
        from pylith.topology.Field import Field
        from pythia.pyre.units.unit import one
        self.vectorFieldType = Field.SCALAR
        self.scale = one / normalizer.getTimeScale()
        self._setComponents(spaceDim)

    def _configure(self):
        """
        Set members based using inventory.
        """
        SolutionSubfield._configure(self)

# FACTORIES ////////////////////////////////////////////////////////////


def soln_subfield():
    """
    Factory associated with SubfieldTraceStrainDot.
    """
    return SubfieldTraceStrainDot()


# End of file