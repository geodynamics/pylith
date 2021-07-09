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
# @file pylith/problems/InitialConditionDomain.py
#
# @brief Python object for specifying initial conditions over the entire domain.
#
# Factory: initial_conditions.

from pylith.problems.InitialCondition import InitialCondition
from .problems import InitialConditionDomain as ModuleInitialCondition


class InitialConditionDomain(InitialCondition, ModuleInitialCondition):
    """Python object for specifying initial conditions over the entire domain.

    FACTORY: initial_conditions
    """

    import pythia.pyre.inventory
    from spatialdata.spatialdb.SimpleDB import SimpleDB

    db = pythia.pyre.inventory.facility("db", family="spatial_database", factory=SimpleDB)
    db.meta["tip"] = "Spatial database with values for initial condition."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="initialconditionsdomain"):
        """Constructor.
        """
        InitialCondition.__init__(self, name)
        return

    def preinitialize(self, problem):
        """Setup initial conditions.
        """
        InitialCondition.preinitialize(self, problem)

        ModuleInitialCondition.setDB(self, self.db)
        return

    def _configure(self):
        """Setup members using inventory.
        """
        InitialCondition._configure(self)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleInitialCondition.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def initial_conditions():
    """Factory associated with InitialConditionDomain.
    """
    return InitialConditionDomain()


# End of file
