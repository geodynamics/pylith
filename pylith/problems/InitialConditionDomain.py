# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

from pylith.problems.InitialCondition import InitialCondition
from .problems import InitialConditionDomain as ModuleInitialCondition


class InitialConditionDomain(InitialCondition, ModuleInitialCondition):
    """
    Initial conditions for the solution over the entire domain.

    Implements `InitialCondition`.
    """
    DOC_CONFIG = {
        "cfg": """
            # Create a single initial condition over the domain.
            [pylithapp.problem]
            ic = [domain]
            ic.domain = pylith.problems.InitialConditionDomain

            [pylithapp.problem.ic.domain]
            db = spatialdata.spatialdb.SimpleGridDB
            db.description = Initial conditions over domain
            db.filename = sheardisp_ic.spatialdb
        """
    }

    import pythia.pyre.inventory
    from spatialdata.spatialdb.SimpleDB import SimpleDB

    db = pythia.pyre.inventory.facility("db", family="spatial_database", factory=SimpleDB)
    db.meta["tip"] = "Spatial database with values for initial condition."

    def __init__(self, name="initialconditionsdomain"):
        """Constructor.
        """
        InitialCondition.__init__(self, name)

    def preinitialize(self, problem):
        """Setup initial conditions.
        """
        InitialCondition.preinitialize(self, problem)

        ModuleInitialCondition.setDB(self, self.db)

    def _configure(self):
        """Setup members using inventory.
        """
        InitialCondition._configure(self)

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
