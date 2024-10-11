# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .ProgressMonitor import ProgressMonitor
from .problems import ProgressMonitorTime as ModuleProgressMonitorTime


class ProgressMonitorTime(ProgressMonitor, ModuleProgressMonitorTime):
    """
    Progress monitor for time-dependent problem.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.timedependent.progress_monitor]
            filename = output/step01-progress.txt
            t_units = year
        """
    }

    import pythia.pyre.inventory

    tUnits = pythia.pyre.inventory.str("t_units", default="year")
    tUnits.meta['tip'] = "Units used for simulation time in output."

    def __init__(self, name="progressmonitortime"):
        """Constructor.
        """
        ProgressMonitor.__init__(self, name)

    def preinitialize(self, defaults):
        """Do minimal initialization.
        """
        ProgressMonitor.preinitialize(self, defaults)
        ModuleProgressMonitorTime.setTimeUnit(self, self.tUnits)

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
