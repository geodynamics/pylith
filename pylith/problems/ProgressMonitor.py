# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.utils.PetscComponent import PetscComponent
from .problems import ProgressMonitor as ModuleProgressMonitor


class ProgressMonitor(PetscComponent, ModuleProgressMonitor):
    """
    Abstract base class for simulation progress monitor.
    """

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="progress.txt")
    filename.meta['tip'] = "Name of output file."

    updatePercent = pythia.pyre.inventory.float(
        "update_percent", default=5.0, validator=pythia.pyre.inventory.greater(0))
    updatePercent.meta['tip'] = "Frequency of progress updates (percent)."

    def __init__(self, name="progressmonitor"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="progress_monitor")

    def preinitialize(self):
        """Do minimal initialization.
        """
        self._createModuleObj()
        ModuleProgressMonitor.setFilename(self, self.filename)
        ModuleProgressMonitor.setUpdatePercent(self, self.updatePercent)
        self._createPath()

    def _createPath(self):
        """Create path for filename if it doesn't exist.
        """
        import os
        import pylith.mpi.mpi as mpi

        relpath = os.path.dirname(self.filename)
        if relpath and not os.path.exists(relpath):
            # Only create directory on master
            isMaster = 0 == mpi.rank()
            if isMaster:
                os.makedirs(relpath)

    def _createModuleObj(self):
        """Create handle to corresponding C++ object.
        """
        raise NotImplementedError("Implement in child class.")


# End of file
