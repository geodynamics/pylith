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


def icFactory(name):
    """
    Factory for initial conditions items.
    """
    from pyre.inventory import facility
    from pylith.problems.InitialConditionDomain import InitialConditionDomain
    return facility(name, family="initial_conditions", factory=InitialConditionDomain)


class TimeDependent(Problem, ModuleTimeDependent):
    """
    Python class for time dependent crustal dynamics problems.

    INVENTORY

    Properties
      - *initial_dt* Initial time step.
      - *start_time* Start time for problem.
      - *end_time* End time for problem.
      - *max_timesteps* Maximum number of time steps.
      - *formulation* Formulation for equation ('implicit' or 'explicit').

    Facilities
      - *initial_conditions* Initial conditions for problem.
      - *progress_monitor* Simple progress monitor via text file.


    FACTORY: problem.
    """

    import pyre.inventory
    from pyre.units.time import year
    from pylith.utils.EmptyBin import EmptyBin

    dtInitial = pyre.inventory.dimensional("initial_dt", default=1.0 * year,
                                           validator=pyre.inventory.greater(0.0 * year))
    dtInitial.meta['tip'] = "Initial time step."

    startTime = pyre.inventory.dimensional("start_time", default=0.0 * year)
    startTime.meta['tip'] = "Start time for problem."

    endTime = pyre.inventory.dimensional("end_time", default=0.1 * year,
                                         validator=pyre.inventory.greaterEqual(0.0 * year))
    endTime.meta['tip'] = "End time for problem."

    maxTimeSteps = pyre.inventory.int("max_timesteps", default=20000, validator=pyre.inventory.greater(0))
    maxTimeSteps.meta['tip'] = "Maximum number of time steps."

    formulation = pyre.inventory.str("formulation", default="implicit",
                                     validator=pyre.inventory.choice(["implicit", "explicit"]))
    formulation.meta['tip'] = "Formulation for equation."

    ic = pyre.inventory.facilityArray("ic", itemFactory=icFactory, factory=EmptyBin)
    ic.meta['tip'] = "Initial conditions."

    shouldNotifyIC = pyre.inventory.bool("notify_observers_ic", default=False)
    shouldNotifyIC.meta["tip"] = "Notify observers of solution with initial conditions."

    from .ProgressMonitorTime import ProgressMonitorTime
    progressMonitor = pyre.inventory.facility(
        "progress_monitor", family="progress_monitor", factory=ProgressMonitorTime)
    progressMonitor.meta['tip'] = "Simple progress monitor via text file."

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

        ModuleTimeDependent.setStartTime(self, self.startTime.value)
        ModuleTimeDependent.setEndTime(self, self.endTime.value)
        ModuleTimeDependent.setInitialTimeStep(self, self.dtInitial.value)
        ModuleTimeDependent.setMaxTimeSteps(self, self.maxTimeSteps)
        if self.formulation == "implicit":
            formulationType = self.IMPLICIT
        elif self.formulation == "explicit":
            formulationType = self.EXPLICIT
        else:
            raise ValueError("Unknown fornulation '{}'.".format(self.formulation))
        ModuleTimeDependent.setFormulation(self, formulationType)
        ModuleTimeDependent.setShouldNotifyIC(self, self.shouldNotifyIC)

        # Preinitialize initial conditions.
        for ic in self.ic.components():
            ic.preinitialize(mesh)
        ModuleTimeDependent.setInitialCondition(self, self.ic.components())

        self.progressMonitor.preinitialize()
        ModuleTimeDependent.setProgressMonitor(self, self.progressMonitor)
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
        if self.startTime > self.endTime:
            raise ValueError("End time {} must be later than start time {}.".format(self.startTime, self.endTime))
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
