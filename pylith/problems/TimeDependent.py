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

from .Problem import Problem
from .problems import TimeDependent as ModuleTimeDependent


def icFactory(name):
    """Factory for initial conditions items.
    """
    from pythia.pyre.inventory import facility
    from pylith.problems.InitialConditionDomain import InitialConditionDomain
    return facility(name, family="initial_conditions", factory=InitialConditionDomain)


class TimeDependent(Problem, ModuleTimeDependent):
    """
    Static, quasistatic, or dynamic time-dependent problem.

    Implements `Problem`.
    """
    DOC_CONFIG = {
        "cfg": """
            # Set boundary conditions, faults, and materials
            bc = [boundary_xpos, boundary_xneg]
            interfaces = [san_andreas, hayward]
            materials = [crust, mantle]

            # Create an initial condition over the domain
            ic = [domain]

            # Turn on gravitational body forces
            gravity_field = spatialdata.spatialdb.GravityField

            # Set the normalizer for nondimensionalizing the problem
            normalizer = spatialdata.units.NondimElasticQuasistatic

            # Set the subfields in the solution
            solution = = pylith.problems.SolnDispLagrange

            # Output the solution for the domain and ground surface
            solution_observers = [domain, ground_surface]

            # Use the quasistatic formulation, linear solver, and set appropriate default solver settings.
            formulation = quasistatic
            solver = linear

            # Use a maximum of 20 time steps to simulation from -0.5 years to 2.0 years with an initial time step of 0.5 years.
            # The first time step will compute the solution at time 0.
            start_time = -0.5*year
            end_time = 2.0*year
            initial_dt = 0.5*year
            max_timesteps = 20

            [pylithapp.greensfns.petsc_defaults]
            solver = True
            monitors = True
        """
    }

    import pythia.pyre.inventory
    from pythia.pyre.units.time import year
    from pylith.utils.EmptyBin import EmptyBin

    dtInitial = pythia.pyre.inventory.dimensional("initial_dt", default=1.0 * year,
                                           validator=pythia.pyre.inventory.greater(0.0 * year))
    dtInitial.meta['tip'] = "Initial time step."

    startTime = pythia.pyre.inventory.dimensional("start_time", default=0.0 * year)
    startTime.meta['tip'] = "Start time for problem."

    endTime = pythia.pyre.inventory.dimensional("end_time", default=0.1 * year,
                                         validator=pythia.pyre.inventory.greaterEqual(0.0 * year))
    endTime.meta['tip'] = "End time for problem."

    maxTimeSteps = pythia.pyre.inventory.int("max_timesteps", default=20000, validator=pythia.pyre.inventory.greater(0))
    maxTimeSteps.meta['tip'] = "Maximum number of time steps."

    ic = pythia.pyre.inventory.facilityArray("ic", itemFactory=icFactory, factory=EmptyBin)
    ic.meta['tip'] = "Initial conditions."

    shouldNotifyIC = pythia.pyre.inventory.bool("notify_observers_ic", default=False)
    shouldNotifyIC.meta["tip"] = "Notify observers of solution with initial conditions."

    from .ProgressMonitorTime import ProgressMonitorTime
    progressMonitor = pythia.pyre.inventory.facility(
        "progress_monitor", family="progress_monitor", factory=ProgressMonitorTime)
    progressMonitor.meta['tip'] = "Simple progress monitor via text file."

    def __init__(self, name="timedependent"):
        """Constructor.
        """
        Problem.__init__(self, name)

    def preinitialize(self, mesh):
        """Setup integrators for each element family (material/quadrature,
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
        ModuleTimeDependent.setShouldNotifyIC(self, self.shouldNotifyIC)

        # Preinitialize initial conditions.
        for ic in self.ic.components():
            ic.preinitialize(mesh)
        ModuleTimeDependent.setInitialCondition(self, self.ic.components())

        self.progressMonitor.preinitialize()
        ModuleTimeDependent.setProgressMonitor(self, self.progressMonitor)

    def run(self, app):
        """Solve time dependent problem.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()

        if 0 == comm.rank:
            self._info.log("Solving problem.")

        ModuleTimeDependent.solve(self)

    def _configure(self):
        """Set members based using inventory.
        """
        Problem._configure(self)
        if self.startTime > self.endTime:
            raise ValueError("End time {} must be later than start time {}.".format(self.startTime, self.endTime))

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleTimeDependent.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def problem():
    """Factory associated with TimeDependent.
    """
    return TimeDependent()


# End of file
