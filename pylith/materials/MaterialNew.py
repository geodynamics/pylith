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

# @file pylith/materials/MaterialNew.py
##
# @brief Python abstract base class for managing physical properties
# and state variables of a material.
##
# Factory: material

from pylith.feassemble.IntegratorPointwise import IntegratorPointwise

# Validator for label


def validateLabel(value):
    """
    Validate descriptive label.
    """
    if 0 == len(value):
        raise ValueError("Descriptive label for material not specified.")
    return value


# MaterialNew class
class MaterialNew(IntegratorPointwise):
    """
    Python material property manager.

    Factory: material
    """

    # INVENTORY //////////////////////////////////////////////////////////

    class Inventory(IntegratorPointwise.Inventory):
        """
        Python object for managing MaterialNew facilities and properties.
        """

        # @class Inventory
        # Python object for managing Material facilities and properties.
        ##
        # \b Properties
        # @li \b id Material identifier (from mesh generator)
        # @li \b label Descriptive label for material.
        ##
        # \b Facilities
        # @li None

        import pyre.inventory

        materialId = pyre.inventory.int("id", default=0)
        materialId.meta['tip'] = "Material identifier (from mesh generator)."

        label = pyre.inventory.str("label", default="", validator=validateLabel)
        label.meta['tip'] = "Descriptive label for material."

        from pylith.meshio.OutputMatElastic import OutputMatElastic
        output = pyre.inventory.facility("output", family="output_manager", factory=OutputMatElastic)
        output.meta['tip'] = "Output manager for elastic material information."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="material"):
        """
        Constructor.
        """
        IntegratorPointwise.__init__(self, name)
        self.output = None
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        try:
            IntegratorPointwise._configure(self)
            self.id(self.inventory.materialId)
            self.label(self.inventory.label)

        except ValueError, err:
            aliases = ", ".join(self.aliases)
            raise ValueError("Error while configuring material (%s):\n%s" % (aliases, err.message))
        return

    def _setupLogging(self):
        """
        Setup event logging.
        """
        if None == self._loggingPrefix:
            self._loggingPrefix = ""

        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.className("FE Material")
        logger.initialize()

        events = ["verify",
                  "init"]
        for event in events:
            logger.registerEvent("%s%s" % (self._loggingPrefix, event))

        self._eventLogger = logger
        return


# End of file
