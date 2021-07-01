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
#

# @file pylith/problems/ProgressMonitorTime.py
##
# @brief Python PyLith object for monitoring progress of time-dependent problem.
##
# Factory: progress_monitor

from .ProgressMonitor import ProgressMonitor
from .problems import ProgressMonitorTime as ModuleProgressMonitorTime


class ProgressMonitorTime(ProgressMonitor, ModuleProgressMonitorTime):
    """Python PyLith object for monitoring progress of time dependent problem.

    Factory: progress_monitor.
    """

    import pythia.pyre.inventory

    tUnits = pythia.pyre.inventory.str("t_units", default="year")
    tUnits.meta['tip'] = "Units for simulation time in output."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="progressmonitortime"):
        """Constructor.
        """
        ProgressMonitor.__init__(self, name)
        return

    def preinitialize(self):
        """Do minimal initialization.
        """
        ProgressMonitor.preinitialize(self)
        ModuleProgressMonitorTime.setTimeUnit(self, self.tUnits)
        return

    # PRIVATE METHODS /////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Create handle to corresponding C++ object.
        """
        ModuleProgressMonitorTime.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def progress_monitor():
    """Factory associated with ProgressMonitorTime.
    """
    return ProgressMonitorTime()


# End of file
