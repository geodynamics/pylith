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
# @file pylith/bc/NeumannTimeDependent.py
#
# @brief Python object for managing a time-dependent Neumann (natural) boundary condition.
#
# Factory: boundary_condition

from .BoundaryCondition import BoundaryCondition
from .bc import NeumannTimeDependent as ModuleNeumannTimeDependent
from pylith.utils.NullComponent import NullComponent


def validateDir(value):
    """Validate direction.
    """
    msg = "Direction must be a 3 component vector (list)."
    if not isinstance(value, list):
        raise ValueError(msg)
    if 3 != len(value):
        raise ValueError(msg)
    try:
        nums = list(map(float, value))
    except:
        raise ValueError(msg)
    return nums


class NeumannTimeDependent(BoundaryCondition, ModuleNeumannTimeDependent):
    """Python object for managing a time-dependent Neumann (natural) boundary condition.

    FACTORY: boundary_condition
    """

    import pythia.pyre.inventory

    scaleName = pythia.pyre.inventory.str("scale_name", default="pressure",
                                   validator=pythia.pyre.inventory.choice(["length", "time", "pressure", "density", "velocity"]))
    scaleName.meta['tip'] = "Type of scale for nondimensionalizing Neumann boundary condition ('pressure' for elasticity)."

    useInitial = pythia.pyre.inventory.bool("use_initial", default=True)
    useInitial.meta['tip'] = "Use initial term in time-dependent expression."

    useRate = pythia.pyre.inventory.bool("use_rate", default=False)
    useRate.meta['tip'] = "Use rate term in time-dependent expression."

    useTimeHistory = pythia.pyre.inventory.bool("use_time_history", default=False)
    useTimeHistory.meta['tip'] = "Use time history term in time-dependent expression."

    dbTimeHistory = pythia.pyre.inventory.facility("time_history", factory=NullComponent, family="temporal_database")
    dbTimeHistory.meta['tip'] = "Time history with normalized amplitude as a function of time."

    refDir1 = pythia.pyre.inventory.list("ref_dir_1", default=[0.0, 0.0, 1.0], validator=validateDir)
    refDir1.meta['tip'] = "First choice for reference direction to discriminate among tangential directions in 3-D."

    refDir2 = pythia.pyre.inventory.list("ref_dir_2", default=[0.0, 1.0, 0.0], validator=validateDir)
    refDir2.meta['tip'] = "Second choice for reference direction to discriminate among tangential directions in 3-D."

    def __init__(self, name="neumanntimedependent"):
        """Constructor.
        """
        BoundaryCondition.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsTimeDependent import AuxSubfieldsTimeDependent
        self.auxiliarySubfields = AuxSubfieldsTimeDependent("auxiliary_subfields")

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log(
                "Performing minimal initialization of time-dependent Neumann boundary condition '%s'." % self.aliases[-1])

        BoundaryCondition.preinitialize(self, problem)

        ModuleNeumannTimeDependent.setScaleName(self, self.scaleName)
        ModuleNeumannTimeDependent.setRefDir1(self, self.refDir1)
        ModuleNeumannTimeDependent.setRefDir2(self, self.refDir2)
        ModuleNeumannTimeDependent.useInitial(self, self.useInitial)
        ModuleNeumannTimeDependent.useRate(self, self.useRate)
        ModuleNeumannTimeDependent.useTimeHistory(self, self.useTimeHistory)
        if not isinstance(self.dbTimeHistory, NullComponent):
            ModuleNeumannTimeDependent.setTimeHistoryDB(self, self.dbTimeHistory)
        return

    def _configure(self):
        """Setup members using inventory.
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
        """Create handle to corresponding C++ object.
        """
        ModuleNeumannTimeDependent.__init__(self)
        return


# Factories

def boundary_condition():
    """Factory associated with NeumannTimeDependent.
    """
    return NeumannTimeDependent()


# End of file
