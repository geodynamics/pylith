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

# @file pylith/problems/ProgressMonitor.py
##
# @brief Python PyLith abstract base class for progress monitor.
##
# Factory: progress_monitor

from pylith.utils.PetscComponent import PetscComponent
from .problems import ProgressMonitor as ModuleProgressMonitor


class ProgressMonitor(PetscComponent, ModuleProgressMonitor):
    """Python abstract base class for progress monitor.

    FACTORY: progress_monitor.
    """

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="progress.txt")
    filename.meta['tip'] = "Name of output file."

    updatePercent = pythia.pyre.inventory.float(
        "update_percent", default=5.0, validator=pythia.pyre.inventory.greater(0))
    updatePercent.meta['tip'] = "Frequency of progress updates (percent)."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="progressmonitor"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="progress_monitor")
        return

    def preinitialize(self):
        """Do minimal initialization.
        """
        self._createModuleObj()
        ModuleProgressMonitor.setFilename(self, self.filename)
        ModuleProgressMonitor.setUpdatePercent(self, self.updatePercent)
        self._createPath()
        return

    # PRIVATE METHODS /////////////////////////////////////////////////////

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
        return

    def _createModuleObj(self):
        """Create handle to corresponding C++ object.
        """
        raise NotImplementedError("Implement in child class.")


# End of file
