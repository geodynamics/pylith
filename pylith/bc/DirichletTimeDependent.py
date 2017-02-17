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
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/bc/DirichletTimeDependent.py
##
# @brief Python object for managing a time-dependent Dirichlet (prescribed
# values) boundary condition.
##
# Factory: boundary_condition

from .DirichletNew import DirichletNew
from .bc import DirichletTimeDependent as ModuleDirichletTimeDependent
from pylith.utils.NullComponent import NullComponent

# DirichletTimeDependent class


class DirichletTimeDependent(DirichletNew,
                             ModuleDirichletTimeDependent):
    """
    Python object for managing a time-dependent Dirichlet (prescribed values)
    boundary condition.

    Factory: boundary_condition
    """

    # INVENTORY //////////////////////////////////////////////////////////

    class Inventory(DirichletNew.Inventory):
        """
        Python object for managing DirichletTimeDependent facilities and properties.
        """

        # @class Inventory
        # Python object for managing ConstraintPointwise facilities and properties.
        ##
        # \b Properties
        # @li \b use_initial Use initial term in time-dependent expression.
        # @li \b use_rate Use rate term in time-dependent expression.
        # @li \b use_time_history Use time history term in time-dependent expression.
        ##
        # \b Facilities
        # @li None

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
        auxFields = pyre.inventory.facilityArray("auxiliary_fields", itemFactory=subfieldFactory, factory=AuxFieldsTimeDependent)
        auxFields.meta['tip'] = "Discretization of constraint parameters."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="dirichlettimedependent"):
        """
        Constructor.
        """
        DirichletNew.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Performing minimal initialization of time-dependent Dirichlet boundary condition '%s'." % self.aliases[-1])

        DirichletNew.preinitialize(self, mesh)

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
        logEvent = "%sverify" % self._loggingPrefix
        self._eventLogger.eventBegin(logEvent)

        DirichletNew.verifyConfiguration(self, self.mesh())
        spaceDim = self.mesh().coordsys().spaceDim()
        for d in self.bcDOF:
            if d < 0 or d >= spaceDim:
                raise ValueError("Attempting to constrain DOF (%d) that doesn't exist for time-dependent Dirichlet boundary condition '%s'. Space dimension is %d." %
                                 (d, self.aliases[-1], spaceDim))

        self._eventLogger.eventEnd(logEvent)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        try:
            if self.inventory.useTimeHistory and isinstance(self.inventory.dbTimeHistory, NullComponent):
                raise ValueError("Missing time history database for time-dependent Dirichlet boundary condition '%s'." % self.aliases[-1])
            if not self.inventory.useTimeHistory and not isinstance(self.inventory.dbTimeHistory, NullComponent):
                self._warning.log("Ignoring time history database setting for time-dependent Dirichlet boundary condition '%s'." % self.aliases[-1])

            DirichletNew._configure(self)
            self.useInitial = self.inventory.useInitial
            self.useRate = self.inventory.useRate
            self.useTimeHistory = self.inventory.useTimeHistory
            self.dbTimeHistory = self.inventory.dbTimeHistory
        except ValueError, err:
            aliases = ", ".join(self.aliases)
            raise ValueError("Error while configuring time-dependent Dirichlet boundary condition "
                             "(%s):\n%s" % (aliases, err.message))
        return

    def _createModuleObj(self):
        """
        Create handle to corresponding C++ object.
        """
        if not hasattr(self, "this"):
            ModuleDirichletTimeDependent.__init__(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def boundary_condition():
    """
    Factory associated with DirichletTimeDependent.
    """
    return DirichletTimeDependent()


# End of file
