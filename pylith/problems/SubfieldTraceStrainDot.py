# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

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