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

# @file pylith/feassemble/ConstraintPointwise.py
##
# @brief Python abstract base class for constraints on operator
# actions with finite-elements.


# ConstraintPointwise class
class ConstraintPointwise(object):
    """
    Python abstract base class for constraints on operator
    actions with finite-elements.
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self):
        """
        Constructor.
        """
        return

    def preinitialize(self, mesh):
        """
        Setup constraint.
        """

        self._setupLogging()
        return

    def verifyConfiguration(self):
        """
        Verify configuration.
        """
        return

    def initialize(self):
        """
        Initialize constraints.
        """
        return

    def finalize(self):
        """
        Cleanup.
        """
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _setupLogging(self):
        """
        Setup event logging.
        """
        if None == self._loggingPrefix:
            self._loggingPrefix = ""

        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.className("FE Constraint")
        logger.initialize()

        events = ["verify",
                  "finalize"]
        for event in events:
            logger.registerEvent("%s%s" % (self._loggingPrefix, event))

        self._eventLogger = logger
        return


# End of file
