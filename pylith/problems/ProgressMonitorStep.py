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
from .problems import ProgressMonitorStep as ModuleProgressMonitorStep


class ProgressMonitorStep(ProgressMonitor, ModuleProgressMonitorStep):
    """
    Progress monitor for problems with a given number of steps, such as Green's functions problem.

    Implementes `ProgressMonitor`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.timedependent.progress_monitor]
            filename = output/greensfns01-progress.txt
        """
    }

    def __init__(self, name="progressmonitorstep"):
        """Constructor.
        """
        ProgressMonitor.__init__(self, name)

    def preinitialize(self, defaults):
        """Do minimal initialization.
        """
        ProgressMonitor.preinitialize(self, defaults)

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
