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
# @file pylith/bc/DirichletTimeDependent.py
#
# @brief Python object for managing a time-dependent Dirichlet (prescribed
# values) boundary condition.
#
# Factory: boundary_condition

from .BoundaryCondition import BoundaryCondition
from .bc import DirichletTimeDependent as ModuleDirichletTimeDependent
from pylith.utils.NullComponent import NullComponent


class DirichletTimeDependent(BoundaryCondition, ModuleDirichletTimeDependent):
    """Python object for managing a time-dependent Dirichlet (prescribed values)
    boundary condition.

    Factory: boundary_condition
    """

    import pythia.pyre.inventory

    constrainedDOF = pythia.pyre.inventory.array(
        "constrained_dof", converter=int, default=[])
    constrainedDOF.meta[
        'tip'] = "Array of constrained degrees of freedom (0=1st DOF, 1=2nd DOF, etc)."

    useInitial = pythia.pyre.inventory.bool("use_initial", default=True)
    useInitial.meta['tip'] = "Use initial term in time-dependent expression."

    useRate = pythia.pyre.inventory.bool("use_rate", default=False)
    useRate.meta['tip'] = "Use rate term in time-dependent expression."

    useTimeHistory = pythia.pyre.inventory.bool(
        "use_time_history", default=False)
    useTimeHistory.meta['tip'] = "Use time history term in time-dependent expression."

    dbTimeHistory = pythia.pyre.inventory.facility(
        "time_history", factory=NullComponent, family="temporal_database")
    dbTimeHistory.meta['tip'] = "Time history with normalized amplitude as a function of time."

    def __init__(self, name="dirichlettimedependent"):
        """Constructor.
        """
        BoundaryCondition.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsTimeDependent import AuxSubfieldsTimeDependent
        self.auxiliarySubfields = AuxSubfieldsTimeDependent(
            "auxiliary_subfields")

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        import numpy

        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log(
                "Performing minimal initialization of time-dependent Dirichlet boundary condition '%s'." % self.aliases[-1])

        BoundaryCondition.preinitialize(self, problem)

        ModuleDirichletTimeDependent.setConstrainedDOF(
            self, numpy.array(self.constrainedDOF, dtype=numpy.int32))
        ModuleDirichletTimeDependent.useInitial(self, self.useInitial)
        ModuleDirichletTimeDependent.useRate(self, self.useRate)
        ModuleDirichletTimeDependent.useTimeHistory(self, self.useTimeHistory)
        if not isinstance(self.dbTimeHistory, NullComponent):
            ModuleDirichletTimeDependent.setTimeHistoryDB(
                self, self.dbTimeHistory)
        return

    def verifyConfiguration(self):
        """Verify compatibility of configuration.
        """
        BoundaryCondition.verifyConfiguration(self, self.mesh())
        spaceDim = self.mesh().coordsys().getSpaceDim()
        for d in self.bcDOF:
            if d < 0 or d >= spaceDim:
                raise ValueError("Attempting to constrain DOF (%d) that doesn't exist for time-dependent Dirichlet boundary condition '%s'. Space dimension is %d." %
                                 (d, self.aliases[-1], spaceDim))
        return

    def _configure(self):
        """Setup members using inventory.
        """
        if 0 == len(self.constrainedDOF):
            raise ValueError("'constrained_dof' must be a zero based integer array of indices corresponding to the "
                             "constrained degrees of freedom.")
        if self.inventory.useTimeHistory and isinstance(self.inventory.dbTimeHistory, NullComponent):
            raise ValueError(
                "Missing time history database for time-dependent Dirichlet boundary condition '%s'." % self.aliases[-1])
        if not self.inventory.useTimeHistory and not isinstance(self.inventory.dbTimeHistory, NullComponent):
            self._warning.log(
                "Ignoring time history database setting for time-dependent Dirichlet boundary condition '%s'." % self.aliases[-1])

        BoundaryCondition._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to corresponding C++ object.
        """
        ModuleDirichletTimeDependent.__init__(self)
        return


# Factories

def boundary_condition():
    """Factory associated with DirichletTimeDependent.
    """
    return DirichletTimeDependent()


# End of file
