# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file pylith/problems/SubfieldPressureDot.py
#
# @brief Python object for pressure_t subfield.
#
# Factory: subfield.

from .SolutionSubfield import SolutionSubfield


class SubfieldPressureDot(SolutionSubfield):
    """
    Object for defining attributes of the time derivative of pressure solution subfield.

    Implements `SolutionSubfield`.
    """
    DOC_CONFIG = {
        "cfg": """
        [pylithapp.problems.solution.subfields.pressure_t]
        alias = pressure_t
        basis_order = 1
        """
    }

    fieldName = "pressure_t"

    def __init__(self, name="subfieldpressure_t"):
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
        self.vectorFieldType = Field.SCALAR
        self.scale = normalizer.getPressureScale() / normalizer.getTimeScale()
        self._setComponents(spaceDim)

    def _configure(self):
        """
        Set members based using inventory.
        """
        SolutionSubfield._configure(self)

# FACTORIES ////////////////////////////////////////////////////////////


def soln_subfield():
    """
    Factory associated with SubfieldPressureDot.
    """
    return SubfieldPressureDot()


# End of file