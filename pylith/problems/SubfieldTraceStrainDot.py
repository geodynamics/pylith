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

    import pythia.pyre.inventory

    from .SolutionSubfield import validateAlias
    userAlias = pythia.pyre.inventory.str(
        "alias", default="trace_strain_t", validator=validateAlias)
    userAlias.meta['tip'] = "Name for subfield."

    fieldName = "trace_strain_t"

    def __init__(self, name="subfieldtracestrain_t"):
        """
        Constructor.
        """
        SolutionSubfield.__init__(self, name)

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