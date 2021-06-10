# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/problems/ProgressMonitorStep.py
##
# @brief Python PyLith object for monitoring progress of Green's functions problem.
##
# Factory: progress_monitor

from .ProgressMonitor import ProgressMonitor
from .problems import ProgressMonitorStep as ModuleProgressMonitorStep


class ProgressMonitorStep(ProgressMonitor, ModuleProgressMonitorStep):
    """Python PyLith object for monitoring progress of Green's functions problem.

    Factory: progress_monitor.
    """

    import pythia.pyre.inventory

    tUnits = pythia.pyre.inventory.str("t_units", default="year")
    tUnits.meta['tip'] = "Units for simulation time in output."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="ProgressMonitorStep"):
        """Constructor.
        """
        ProgressMonitor.__init__(self, name)
        return

    def preinitialize(self):
        """Do minimal initialization.
        """
        ProgressMonitor.preinitialize(self)
        ModuleProgressMonitorStep.setTimeUnit(self, self.tUnits)
        return

    # PRIVATE METHODS /////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Create handle to corresponding C++ object.
        """
        ModuleProgressMonitorStep.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def progress_monitor():
    """Factory associated with ProgressMonitorStep.
    """
    return ProgressMonitorStep()


# End of file
