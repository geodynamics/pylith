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

# @file pylith/problems/SubfieldTemperature.py
##
# @brief Python object for temperature subfield.
##
# Factory: subfield.

from .SolutionSubfield import SolutionSubfield

# SubfieldTemperature class


class SubfieldTemperature(SolutionSubfield):
    """
    Python object for temperature subfield.

    Factory: subfield.
    """

    # INVENTORY //////////////////////////////////////////////////////////

    class Inventory(SolutionSubfield.Inventory):
        """
        Python object for managing SubfieldTemperature facilities and properties.
        """

        # @class Inventory
        # Python object for managing SubfieldTemperature facilities and properties.
        ##
        # \b Properties
        # @li \b name Name for subfield.
        ##
        # \b Facilities
        # @li None

        import pyre.inventory

        from .SolutionSubfield import validateName
        fieldName = pyre.inventory.str("name", default="temperature", validator=validateName)
        fieldName.meta['tip'] = "Name for subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="subfieldtemperature"):
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
        self.vectorFieldType = Field.SCALAR
        self.scale = normalizer.temperatureScale()
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
    Factory associated with SubfieldTemperature.
    """
    return SubfieldTemperature()


# End of file
