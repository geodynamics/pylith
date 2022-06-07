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