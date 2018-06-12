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
# @file pylith/problems/TimeDependent.py
#
# @brief Python class for time dependent crustal
# dynamics problems.
#
# Factory: problem.

from .Problem import Problem
from .problems import TimeDependent as ModuleTimeDependent


class TimeDependent(Problem, ModuleTimeDependent):
    """
    Python class for time dependent crustal dynamics problems.

    INVENTORY

    Properties
      - *initial_dt* Initial time step.
      - *start_time* Start time for problem.
      - *total_time* Time duration of problem.
      - *max_timesteps* Maximum number of time steps.

    Facilities
      - *initializer* Problem initializer.
      - *progress_monitor* Simple progress monitor via text file.


    FACTORY: problem.
    """

    # INVENTORY //////////////////////////////////////////////////////////
    #

    import pyre.inventory
    from pyre.units.time import year

    dtInitial = pyre.inventory.dimensional("initial_dt", default=1.0 * year, validator=pyre.inventory.greater(0.0 * year))
    dtInitial.meta['tip'] = "Initial time step."

    startTime = pyre.inventory.dimensional("start_time", default=0.0 * year)
    startTime.meta['tip'] = "Start time for problem."

    totalTime = pyre.inventory.dimensional("total_time", default=0.1 * year, validator=pyre.inventory.greaterEqual(0.0 * year))
    totalTime.meta['tip'] = "Time duration of problem."

    maxTimeSteps = pyre.inventory.int("max_timesteps", default=20000, validator=pyre.inventory.greater(0))
    maxTimeSteps.meta['tip'] = "Maximum number of time steps."

    #from ProgressMonitorTime import ProgressMonitorTime
    #progressMonitor = pyre.inventory.facility("progress_monitor", family="progress_monitor", factory=ProgressMonitorTime)
    #progressMonitor.meta['tip'] = "Simple progress monitor via text file."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="timedependent"):
        """
        Constructor.
        """
        Problem.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Setup integrators for each element family (material/quadrature,
        bc/quadrature, etc.).
        """
        self._setupLogging()

        import weakref
        self.mesh = weakref.ref(mesh)

        Problem.preinitialize(self, mesh)

        ModuleTimeDependent.startTime(self, self.startTime.value)
        ModuleTimeDependent.dtInitial(self, self.dtInitial.value)
        ModuleTimeDependent.totalTime(self, self.totalTime.value)
        ModuleTimeDependent.maxTimeSteps(self, self.maxTimeSteps)

        return

    def run(self, app):
        """
        Solve time dependent problem.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()

        if 0 == comm.rank:
            self._info.log("Solving problem.")

        ModuleTimeDependent.solve(self)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        Problem._configure(self)
        return

    def _createModuleObj(self):
        """
        Create handle to C++ object.
        """
        ModuleTimeDependent.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def problem():
    """
    Factory associated with TimeDependent.
    """
    return TimeDependent()


# End of file
