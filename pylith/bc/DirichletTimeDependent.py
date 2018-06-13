# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/bc/DirichletTimeDependent.py
#
# @brief Python object for managing a time-dependent Dirichlet (prescribed
# values) boundary condition.
#
# Factory: boundary_condition

from .Dirichlet import Dirichlet
from .bc import DirichletTimeDependent as ModuleDirichletTimeDependent
from pylith.utils.NullComponent import NullComponent


class DirichletTimeDependent(Dirichlet,
                             ModuleDirichletTimeDependent):
    """
    Python object for managing a time-dependent Dirichlet (prescribed values)
    boundary condition.

    INVENTORY

    Properties
      - *use_initial* Use initial term in time-dependent expression.
      - *use_rate* Use rate term in time-dependent expression.
      - *use_time_history* Use time history term in time-dependent expression.

    Facilities
      - None

    Factory: boundary_condition
    """

    import pyre.inventory

    useInitial = pyre.inventory.bool("use_initial", default=True)
    useInitial.meta['tip'] = "Use initial term in time-dependent expression."

    useRate = pyre.inventory.bool("use_rate", default=False)
    useRate.meta['tip'] = "Use rate term in time-dependent expression."

    useTimeHistory = pyre.inventory.bool("use_time_history", default=False)
    useTimeHistory.meta['tip'] = "Use time history term in time-dependent expression."

    dbTimeHistory = pyre.inventory.facility("time_history", factory=NullComponent, family="temporal_database")
    dbTimeHistory.meta['tip'] = "Time history with normalized amplitude as a function of time."

    from .AuxFieldsTimeDependent import AuxFieldsTimeDependent
    from pylith.topology.AuxSubfield import subfieldFactory
    auxSubfields = pyre.inventory.facilityArray("auxiliary_subfields", itemFactory=subfieldFactory, factory=AuxFieldsTimeDependent)
    auxSubfields.meta['tip'] = "Discretization of constraint parameters."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="dirichlettimedependent"):
        """
        Constructor.
        """
        Dirichlet.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Performing minimal initialization of time-dependent Dirichlet boundary condition '%s'." % self.aliases[-1])

        Dirichlet.preinitialize(self, mesh)

        ModuleDirichletTimeDependent.useInitial(self, self.useInitial)
        ModuleDirichletTimeDependent.useRate(self, self.useRate)
        ModuleDirichletTimeDependent.useTimeHistory(self, self.useTimeHistory)
        if not isinstance(self.dbTimeHistory, NullComponent):
            ModuleDirichletTimeDependent.dbTimeHistory(self.dbTimeHistory)
        return

    def verifyConfiguration(self):
        """
        Verify compatibility of configuration.
        """
        Dirichlet.verifyConfiguration(self, self.mesh())
        spaceDim = self.mesh().coordsys().spaceDim()
        for d in self.bcDOF:
            if d < 0 or d >= spaceDim:
                raise ValueError("Attempting to constrain DOF (%d) that doesn't exist for time-dependent Dirichlet boundary condition '%s'. Space dimension is %d." %
                                 (d, self.aliases[-1], spaceDim))
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        if self.inventory.useTimeHistory and isinstance(self.inventory.dbTimeHistory, NullComponent):
            raise ValueError("Missing time history database for time-dependent Dirichlet boundary condition '%s'." % self.aliases[-1])
        if not self.inventory.useTimeHistory and not isinstance(self.inventory.dbTimeHistory, NullComponent):
            self._warning.log("Ignoring time history database setting for time-dependent Dirichlet boundary condition '%s'." % self.aliases[-1])

        Dirichlet._configure(self)
        return

    def _createModuleObj(self):
        """
        Create handle to corresponding C++ object.
        """
        ModuleDirichletTimeDependent.__init__(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def boundary_condition():
    """
    Factory associated with DirichletTimeDependent.
    """
    return DirichletTimeDependent()


# End of file
