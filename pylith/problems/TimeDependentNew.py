#!/usr/bin/env python
#
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

# @file pylith/problems/TimeDependentNew.py
##
# @brief Python class for time dependent crustal
# dynamics problems.
##
# Factory: problem.

from ProblemNew import ProblemNew
from problems import TimeDependent as ModuleTimeDependent

# TimeDependentNew class


class TimeDependentNew(ProblemNew, ModuleTimeDependent):
    """
    Python class for time dependent crustal dynamics problems.

    Factory: problem.
    """

    # INVENTORY //////////////////////////////////////////////////////////

    class Inventory(ProblemNew.Inventory):
        """
        Python object for managing TimeDependentNew facilities and properties.
        """

        # @class Inventory
        # Python object for managing TimeDependentNew facilities and properties.
        ##
        # \b Properties
        # @li \b initial_dt Initial time step.
        # @li \b start_time Start time for problem.
        # @li \b total_time Time duration of problem.
        # @li \b max_timesteps Maximum number of time steps.
        ##
        # \b Facilities
        # @li \b initializer Problem initializer.
        # @li \b progress_monitor Simple progress monitor via text file.

        import pyre.inventory
        from pyre.units.time import year

        dtInitial = pyre.inventory.dimensional("initial_dt", default=1.0 * year, validator=pyre.inventory.greater(0.0 * year))
        dtInitial.meta['tip'] = "Initial time step."

        startTime = pyre.inventory.dimensional("start_time", default=0.0 * year)
        startTime.meta['tip'] = "Start time for problem."

        totalTime = pyre.inventory.dimensional("total_time", default=0.0 * year, validator=pyre.inventory.greaterEqual(0.0 * year))
        totalTime.meta['tip'] = "Time duration of problem."

        maxTimeSteps = pyre.inventory.int("max_timesteps", default=20000, validator=pyre.inventory.greater(0))
        maxTimeSteps.meta['tip'] = "Maximum number of time steps."

        from ProgressMonitorTime import ProgressMonitorTime
        progressMonitor = pyre.inventory.facility("progress_monitor", family="progress_monitor", factory=ProgressMonitorTime)
        progressMonitor.meta['tip'] = "Simple progress monitor via text file."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="timedependent"):
        """
        Constructor.
        """
        ProblemNew.__init__(self, name)
        ModuleTimeDependent.__init__(self)
        self._loggingPrefix = "PrTD "
        return

    def preinitialize(self, mesh):
        """
        Setup integrators for each element family (material/quadrature,
        bc/quadrature, etc.).
        """
        self._setupLogging()
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()

        if 0 == comm.rank:
            self._info.log("Pre-initializing problem.")
        import weakref
        self.mesh = weakref.ref(mesh)
        ProblemNew.preinitialize(self, mesh)
        return

    def verifyConfiguration(self):
        """
        Verify compatibility of configuration.
        """
        ProblemNew.verifyConfiguration(self)
        return

    def initialize(self):
        """
        Setup integrators for each element family (material/quadrature,
        bc/quadrature, etc.).
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()

        if 0 == comm.rank:
            self._info.log("Initializing problem.")
        ProblemNew.initialize(self)
        return

    def run(self, app):
        """
        Solve time dependent problem.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()

        if 0 == comm.rank:
            self._info.log("Solving problem.")
        self.checkpointTimer.toplevel = app  # Set handle for saving state

        ModuleTimeDependent.create(self)
        ModuleTimeDependent.initialize()
        ModuleTimeDependent.solve()
        return

    def finalize(self):
        """
        Cleanup after running problem.
        """
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        ProblemNew._configure(self)
        ModuleTimeDependent.startTime(self, self.inventory.startTime.value)
        ModuleTimeDependent.dtInitial(self, self.inventory.dtInitial.value)
        ModuleTimeDependent.totalTime(self, self.inventory.totalTime.value)
        ModuleTimeDependent.maxTimeSteps(self, self.inventory.maxTimeSteps)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def problem():
    """
    Factory associated with TimeDependentNew.
    """
    return TimeDependentNew()


# End of file
