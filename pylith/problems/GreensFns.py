# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/problems/GreensFns.py
#
# @brief Python class for Green's functions problem
#
# Factory: problem.

from .Problem import Problem
from .problems import GreensFns as ModuleGreensFns


class GreensFns(Problem, ModuleGreensFns):
    """Python class for Green's functions problem.

    FACTORY: problem.
    """

    import pythia.pyre.inventory

    faultId = pythia.pyre.inventory.int("fault_id", default=100)
    faultId.meta['tip'] = "Id of fault on which to impose impulses."

    from .ProgressMonitorStep import ProgressMonitorStep
    progressMonitor = pythia.pyre.inventory.facility(
        "progress_monitor", family="progress_monitor", factory=ProgressMonitorStep)
    progressMonitor.meta['tip'] = "Simple progress monitor via text file."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="greensfns"):
        """Constructor.
        """
        Problem.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """Setup integrators for each element family (material/quadrature,
        bc/quadrature, etc.).
        """
        self._setupLogging()

        import weakref
        self.mesh = weakref.ref(mesh)

        Problem.preinitialize(self, mesh)

        ModuleGreensFns.setFaultId(self, self.faultId)

        self.progressMonitor.preinitialize()
        ModuleGreensFns.setProgressMonitor(self, self.progressMonitor)
        return

    def run(self, app):
        """Solve time dependent problem.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()

        if 0 == comm.rank:
            self._info.log("Solving problem.")

        ModuleGreensFns.solve(self)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        Problem._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleGreensFns.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def problem():
    """Factory associated with GreensFns.
    """
    return GreensFns()


# End of file
