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

    def preinitialize(self):
        """Do minimal initialization.
        """
        ProgressMonitor.preinitialize(self)

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
