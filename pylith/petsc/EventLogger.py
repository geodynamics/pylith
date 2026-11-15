# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from .petsc import EventLogger as ModuleEventLogger


class EventLogger(ModuleEventLogger):
    """
    Manage event logging using PETSc.

    Each logger object manages the events for a single 'logging class'.
    """

    def __init__(self):
        """Constructor."""
        ModuleEventLogger.__init__(self)
        self.events = {}  # dict of events with counts for current logging.

    def stagePush(self, name):
        """Log stage begin."""
        ModuleEventLogger.stagePush(self, self.getStageId(name))


# End of file
