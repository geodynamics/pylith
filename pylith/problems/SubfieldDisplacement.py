#!/usr/bin/env python
#
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

# @file pylith/problems/SubfieldDisplacement.py
##
# @brief Python object for displacement subfield.
##
# Factory: subfield.

from .SolutionSubfield import SolutionSubfield

# SubfieldDisplacement class


class SubfieldDisplacement(SolutionSubfield):
    """
    Python object for displacement subfield.

    Factory: subfield.
    """

    # INVENTORY //////////////////////////////////////////////////////////

    class Inventory(SolutionSubfield.Inventory):
        """
        Python object for managing SubfieldDisplacement facilities and properties.
        """

        # @class Inventory
        # Python object for managing SubfieldDisplacement facilities and properties.
        ##
        # \b Properties
        # @li \b name Name for subfield.
        ##
        # \b Facilities
        # @li None

        import pyre.inventory

        from .SolutionSubfield import validateName
        fieldName = pyre.inventory.str("name", default="displacement", validator=validateName)
        fieldName.meta['tip'] = "Name for subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="subfielddisplacement"):
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
        self.vectorFieldType = Field.VECTOR
        self.scale = normalizer.lengthScale()
        self._setComponents(spaceDim)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def soln_subfield():
    """
    Factory associated with SubfieldDisplacement.
    """
    return SubfieldDisplacement()


# End of file
