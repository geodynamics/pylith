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
# @file pylith/bc/NeumannTimeDependent.py
#
# @brief Python object for managing a time-dependent Neumann (natural) boundary condition.
#
# Factory: boundary_condition

from .Neumann import Neumann
from .bc import NeumannTimeDependent as ModuleNeumannTimeDependent
from pylith.utils.NullComponent import NullComponent

# NeumannTimeDependent class


class NeumannTimeDependent(Neumann,
                           ModuleNeumannTimeDependent):
    """
    Python object for managing a time-dependent Neumann (prescribed values)
    boundary condition.

    Factory: boundary_condition
    """

    # INVENTORY //////////////////////////////////////////////////////////
    #
    # \b Properties
    # @li \b use_initial Use initial term in time-dependent expression.
    # @li \b use_rate Use rate term in time-dependent expression.
    # @li \b use_time_history Use time history term in time-dependent expression.
    #
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
    auxSubfields = pyre.inventory.facilityArray("auxiliary_subfields", itemFactory=subfieldFactory, factory=AuxFieldsTimeDependent)
    auxSubfields.meta['tip'] = "Discretization of constraint parameters."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="neumanntimedependent"):
        """
        Constructor.
        """
        Neumann.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Performing minimal initialization of time-dependent Neumann boundary condition '%s'." % self.aliases[-1])

        Neumann.preinitialize(self, mesh)

        ModuleNeumannTimeDependent.useInitial(self, self.useInitial)
        ModuleNeumannTimeDependent.useRate(self, self.useRate)
        ModuleNeumannTimeDependent.useTimeHistory(self, self.useTimeHistory)
        if not isinstance(self.dbTimeHistory, NullComponent):
            ModuleNeumannTimeDependent.dbTimeHistory(self.dbTimeHistory)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        if self.inventory.useTimeHistory and isinstance(self.inventory.dbTimeHistory, NullComponent):
            raise ValueError("Missing time history database for time-dependent Neumann boundary condition '%s'." % self.aliases[-1])
        if not self.inventory.useTimeHistory and not isinstance(self.inventory.dbTimeHistory, NullComponent):
            self._warning.log("Ignoring time history database setting for time-dependent Neumann boundary condition '%s'." % self.aliases[-1])

        Neumann._configure(self)
        return

    def _createModuleObj(self):
        """
        Create handle to corresponding C++ object.
        """
        ModuleNeumannTimeDependent.__init__(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def boundary_condition():
    """
    Factory associated with NeumannTimeDependent.
    """
    return NeumannTimeDependent()


# End of file
