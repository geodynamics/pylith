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

from .utils import EventLogger as ModuleEventLogger


class EventLogger(ModuleEventLogger):
    """
    Manage event logging using PETSc.

    Each logger object manages the events for a single 'logging class'.
    """

    def __init__(self):
        """Constructor.
        """
        ModuleEventLogger.__init__(self)
        self.events = {}  # dict of events with counts for current logging.

    def registerEvent(self, name):
        """Register event.
        """
        self.events[name] = 0  # Set log count to 0
        return ModuleEventLogger.registerEvent(self, name)

    def eventBegin(self, name):
        """Log event begin.
        """
        if self.events[name] == 0:  # prevent double counting
            ModuleEventLogger.eventBegin(self, self.getEventId(name))
        self.events[name] += 1

    def eventEnd(self, name):
        """Log event end.
        """
        if self.events[name] > 0:
            self.events[name] -= 1
        if 0 == self.events[name]:  # prevent double counting
            ModuleEventLogger.eventEnd(self, self.getEventId(name))

    def stagePush(self, name):
        """Log stage begin.
        """
        ModuleEventLogger.stagePush(self, self.getStageId(name))


# End of file
