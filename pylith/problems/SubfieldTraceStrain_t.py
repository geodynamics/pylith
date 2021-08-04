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
#

# @file pylith/problems/SubfieldTraceStrain_t.py
##
# @brief Python object for trace_strain_t subfield.
##
# Factory: subfield.

from .SolutionSubfield import SolutionSubfield


class SubfieldTraceStrain_t(SolutionSubfield):
    """
    Python object for trace_strain_t subfield.

    INVENTORY

    Properties
      - *alias* User-specified name for subfield.

    Facilities
      - None

    FACTORY: subfield
    """

    import pythia.pyre.inventory

    from .SolutionSubfield import validateAlias
    userAlias = pythia.pyre.inventory.str(
        "alias", default="trace_strain_t", validator=validateAlias)
    userAlias.meta['tip'] = "Name for subfield."

    fieldName = "trace_strain_t"

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="subfieldtracestrain_t"):
        """
        Constructor.
        """
        SolutionSubfield.__init__(self, name)
        return

    def initialize(self, normalizer, spaceDim):
        """
        Initialize subfield metadata.
        """
        from pylith.topology.Field import Field
        from pythia.pyre.units.unit import one
        self.vectorFieldType = Field.SCALAR
        self.scale = one / normalizer.getTimeScale()
        self._setComponents(spaceDim)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        SolutionSubfield._configure(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def soln_subfield():
    """
    Factory associated with SubfieldTraceStrain_t.
    """
    return SubfieldTraceStrain_t()


# End of file
