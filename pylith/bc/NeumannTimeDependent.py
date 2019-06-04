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

from .BoundaryCondition import BoundaryCondition
from .bc import NeumannTimeDependent as ModuleNeumannTimeDependent
from pylith.utils.NullComponent import NullComponent


class NeumannTimeDependent(BoundaryCondition, ModuleNeumannTimeDependent):
    """
    Python object for managing a time-dependent Neumann (natural) boundary condition.

    INVENTORY

    Properties
      - *scale_name* Type of scale for nondimenaionlizing Neumann boundary condition (e.g., "pressure" for elasticity").
      - *use_initial* Use initial term in time-dependent expression.
      - *use_rate* Use rate term in time-dependent expression.
      - *use_time_history* Use time history term in time-dependent expression.

    Facilities
      - *auxiliary_subfields* Discretization of time-dependent Neumann parameters.

    FACTORY: boundary_condition
    """

    import pyre.inventory

    scaleName = pyre.inventory.str("scale_name", default="pressure",
                                   validator=pyre.inventory.choice(["length", "time", "pressure", "density", "velocity"]))
    scaleName.meta['tip'] = "Type of scale for nondimensionalizing Neumann boundary condition ('pressure' for elasticity)."

    useInitial = pyre.inventory.bool("use_initial", default=True)
    useInitial.meta['tip'] = "Use initial term in time-dependent expression."

    useRate = pyre.inventory.bool("use_rate", default=False)
    useRate.meta['tip'] = "Use rate term in time-dependent expression."

    useTimeHistory = pyre.inventory.bool("use_time_history", default=False)
    useTimeHistory.meta['tip'] = "Use time history term in time-dependent expression."

    dbTimeHistory = pyre.inventory.facility("time_history", factory=NullComponent, family="temporal_database")
    dbTimeHistory.meta['tip'] = "Time history with normalized amplitude as a function of time."

    from .AuxSubfieldsTimeDependent import AuxSubfieldsTimeDependent
    from pylith.topology.Subfield import subfieldFactory
    auxiliarySubfields = pyre.inventory.facilityArray(
        "auxiliary_subfields", itemFactory=subfieldFactory, factory=AuxSubfieldsTimeDependent)
    auxiliarySubfields.meta['tip'] = "Discretization of time-dependent Neumann parameters."

    def __init__(self, name="neumanntimedependent"):
        """
        Constructor.
        """
        BoundaryCondition.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log(
                "Performing minimal initialization of time-dependent Neumann boundary condition '%s'." % self.aliases[-1])

        BoundaryCondition.preinitialize(self, mesh)

        ModuleNeumannTimeDependent.setScaleName(self, self.scaleName)
        ModuleNeumannTimeDependent.useInitial(self, self.useInitial)
        ModuleNeumannTimeDependent.useRate(self, self.useRate)
        ModuleNeumannTimeDependent.useTimeHistory(self, self.useTimeHistory)
        if not isinstance(self.dbTimeHistory, NullComponent):
            ModuleNeumannTimeDependent.dbTimeHistory(self.dbTimeHistory)
        return

    def _configure(self):
        """
        Setup members using inventory.
        """
        if self.inventory.useTimeHistory and isinstance(self.inventory.dbTimeHistory, NullComponent):
            raise ValueError(
                "Missing time history database for time-dependent Neumann boundary condition '%s'." % self.aliases[-1])
        if not self.inventory.useTimeHistory and not isinstance(self.inventory.dbTimeHistory, NullComponent):
            self._warning.log(
                "Ignoring time history database setting for time-dependent Neumann boundary condition '%s'." % self.aliases[-1])

        BoundaryCondition._configure(self)
        return

    def _createModuleObj(self):
        """
        Create handle to corresponding C++ object.
        """
        ModuleNeumannTimeDependent.__init__(self)
        return


# Factories

def boundary_condition():
    """
    Factory associated with NeumannTimeDependent.
    """
    return NeumannTimeDependent()


# End of file
