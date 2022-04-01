# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

from .SolutionSubfield import SolutionSubfield


class SubfieldTemperature(SolutionSubfield):
    """
    Object for defining attributes of the temperature solution subfield.

    Implements `SolutionSubfield`.
    """
    DOC_CONFIG = {
        "cfg": """
        [pylithapp.problems.solution.subfields.temperature]
        alias = absolute_temperature
        basis_order = 1
        """
    }

    import pythia.pyre.inventory

    from .SolutionSubfield import validateAlias
    userAlias = pythia.pyre.inventory.str("alias", default="temperature", validator=validateAlias)
    userAlias.meta['tip'] = "Name for subfield."

    fieldName = "temperature"

    def __init__(self, name="subfieldtemperature"):
        """Constructor.
        """
        SolutionSubfield.__init__(self, name)

    def initialize(self, normalizer, spaceDim):
        """Initialize subfield metadata.
        """
        from pylith.topology.Field import Field
        self.vectorFieldType = Field.SCALAR
        self.scale = normalizer.getTemperatureScale()
        self._setComponents(spaceDim)

    def _configure(self):
        """Set members based using inventory.
        """
        SolutionSubfield._configure(self)

# FACTORIES ////////////////////////////////////////////////////////////


def soln_subfield():
    """Factory associated with SubfieldTemperature.
    """
    return SubfieldTemperature()


# End of file
