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
    """
    Static Green's function problem type with each Green's function corresponding to a fault slip impulses.

    Implements `Problem`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp]
            problem = pylith.problems.GreensFns

            [pylithapp.greensfns]
            label = fault
            label_value = 1

            # Set appropriate default solver settings.
            set_solver_defaults = True

            interfaces = [fault]
            interfaces.fault = pylith.faults.FaultCohesiveImpulses

            [pylithapp.greensfns.interfaces.fault]
            label = fault
            label_value = 20
            
            # Impulses for left-lateral slip (dof=1)
            impulse_dof = [1]
            threshold = 0.5

            # Create impulses at all points on the fault by specifying a uniform amplitude of 1.0.
            # Impulses will be applied at any location with a slip component greater than the threshold.
            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Slip impulse amplitude
            db_auxiliary_field.values = [slip_left_lateral, slip_opening]
            db_auxiliary_field.data = [1.0*m, 0.0*m]

            # Represent the impulse as a linear variation in slip centered on each point.
            auxiliary_subfields.slip.basis_order = 1

            [pylithapp.greensfns.petsc_defaults]
            solver = True
            monitors = True
        """
    }

    import pythia.pyre.inventory

    faultLabelName = pythia.pyre.inventory.str("label", default="fault")
    faultLabelName.meta['tip'] = "Name of label identifier for fault surface on which to impose impulses."

    faultLabelValue = pythia.pyre.inventory.int("label_value", default=1)
    faultLabelValue.meta['tip'] = "Value of label identifier for fault surface on which to impose impulses."

    from .ProgressMonitorStep import ProgressMonitorStep
    progressMonitor = pythia.pyre.inventory.facility(
        "progress_monitor", family="progress_monitor", factory=ProgressMonitorStep)
    progressMonitor.meta['tip'] = "Simple progress monitor via text file."

    def __init__(self, name="greensfns"):
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

        ModuleGreensFns.setFaultLabelName(self, self.faultLabelName)
        ModuleGreensFns.setFaultLabelValue(self, self.faultLabelValue)

        self.progressMonitor.preinitialize()
        ModuleGreensFns.setProgressMonitor(self, self.progressMonitor)

    def run(self, app):
        """Solve time dependent problem.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()

        if 0 == comm.rank:
            self._info.log("Solving problem.")

        ModuleGreensFns.solve(self)

    def _configure(self):
        """Set members based using inventory.
        """
        Problem._configure(self)

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleGreensFns.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def problem():
    """Factory associated with GreensFns.
    """
    return GreensFns()


# End of file
