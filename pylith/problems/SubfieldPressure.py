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
from spatialdata.units.ElasticityScales import ElasticityScales


class SubfieldPressure(SolutionSubfield):
    """
    Object for defining attributes of the pressure solution subfield.

    Implements `SolutionSubfield`.
    """

    DOC_CONFIG = {
        "cfg": """
        [pylithapp.problems.solution.subfields.pressure]
        alias = pressure
        basis_order = 1
        """
    }

    fieldName = "pressure"

    def __init__(self, name="subfieldpressure"):
        """Constructor."""
        SolutionSubfield.__init__(self, name)
        return

    def _defaults(self):
        self.userAlias = self.fieldName

    def initialize(self, scales, spaceDim):
        """Initialize subfield metadata."""
        from pylith.topology.Field import Field

        self.vectorFieldType = Field.SCALAR
        self.scale = ElasticityScales.getFluidPressureScale(scales)
        self._setComponents(spaceDim)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory."""
        SolutionSubfield._configure(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////


def soln_subfield():
    """Factory associated with SubfieldPressure."""
    return SubfieldPressure()


# End of file
