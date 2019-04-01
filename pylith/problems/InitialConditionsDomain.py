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
# @file pylith/problems/InitialConditionsDomain.py
#
# @brief Python object for specifying initial conditions over the entire domain.
#
# Factory: initial_conditions.

from pylith.problems.InitialConditions import InitialConditions
from .problems import InitialConditionsDomain as ModuleInitialConditions


class InitialConditionsDomain(InitialConditions, ModuleInitialConditions):
    """
    Python object for specifying initial conditions over the entire domain.

    INVENTORY

    Properties
      - None

    Facilities
      - *db* Spatial database with values for initial conditions.
    """

    import pyre.inventory

    db = pyre.inventory.facility("db", default=SimpleDB, factory="spatial_database")
    db.meta["tip"] = "Spatial database with values for initial condition."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="initialconditionsdomain"):
        """
        Constructor.
        """
        InitialConditions.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Setup initial conditions.
        """
        InitialConditions.preinitialize(self, mesh)

        ModuleInitialConditions.setDB(self, self.db)
        return

    def _configure(self):
        """
        Setup members using inventory.
        """
        InitialConditions._configure(self)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        ModuleInitialConditionsPatch.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def initial_conditions():
    """
    Factory associated with InitialConditionsDomain.
    """
    return InitialConditionsDomain()


# End of file
