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


def icFactory(name):
    """Factory for initial conditions items.
    """
    from pythia.pyre.inventory import facility
    from pylith.problems.InitialConditionDomain import InitialConditionDomain
    return facility(name, family="initial_conditions", factory=InitialConditionDomain)


class GreensFns(Problem, ModuleGreensFns):
    """Python class for Green's functions problem.

    FACTORY: problem.
    """

    import pythia.pyre.inventory
    from pythia.pyre.units.time import year
    from pylith.utils.EmptyBin import EmptyBin

    faultId = pythia.pyre.inventory.dimensional("fault_id", default=100)
    faultId.meta['tip'] = "Id of fault on which to impose impulses."

    numImpulses = pythia.pyre.inventory.dimensional("num_impulses", default=1, validator=pythia.pyre.inventory.greater(0))
    numImpulses.meta['tip'] = "Number of impulses."

    ic = pythia.pyre.inventory.facilityArray("ic", itemFactory=icFactory, factory=EmptyBin)
    ic.meta['tip'] = "Initial conditions."

    shouldNotifyIC = pythia.pyre.inventory.bool("notify_observers_ic", default=False)
    shouldNotifyIC.meta["tip"] = "Notify observers of solution with initial conditions."

    from .ProgressMonitorTime import ProgressMonitorTime
    progressMonitor = pythia.pyre.inventory.facility(
        "progress_monitor", family="progress_monitor", factory=ProgressMonitorTime)
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

        ModuleGreensFns.setfaultId(self, self.faultId.value)
        ModuleGreensFns.setnumImpulses(self, self.numImpulses.value)
        ModuleGreensFns.setShouldNotifyIC(self, self.shouldNotifyIC)

        # Preinitialize initial conditions.
        for ic in self.ic.components():
            ic.preinitialize(mesh)
        ModuleGreensFns.setInitialCondition(self, self.ic.components())

        # Find fault for impulses
        found = False
        for fault in self.interfaces.components():
            if self.faultId == fault.id():
                self.source = fault
                found = True
                break
        if not found:
            raise ValueError("Could not find fault interface with id '%d' for "
                            "Green's function impulses." % self.faultId)

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
