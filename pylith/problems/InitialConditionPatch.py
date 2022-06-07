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
from .problems import InitialConditionPatch as ModuleInitialCondition


class InitialConditionPatch(InitialCondition, ModuleInitialCondition):
    """
    Initial conditions over a portion of the domain (patch).

    Implements `InitialCondition`.
    """
    DOC_CONFIG = {
        "cfg": """
            # Create separate initial conditions for two materials.
            # This is often useful if the materials have different properties.
            [pylithapp.problem]
            ic = [mat1, mat2]
            ic.mat1 = pylith.problems.InitialConditionPatch
            ic.mat2 = pylith.problems.InitialConditionPatch

            [pylithapp.problem.ic.mat1]
            label_value = 1
            db = spatialdata.spatialdb.SimpleGridDB
            db.description = Initial conditions over material 1
            db.filename = shearmat1_ic.spatialdb

            [pylithapp.problem.ic.mat2]
            label_value = 2
            db = spatialdata.spatialdb.SimpleGridDB
            db.description = Initial conditions over material 2
            db.filename = shearmat2_ic.spatialdb
        """
    }

    import pythia.pyre.inventory
    from spatialdata.spatialdb.SimpleDB import SimpleDB

    labelName = pythia.pyre.inventory.str("label", default="material-id")
    labelName.meta['tip'] = "Name of label for patch."

    labelValue = pythia.pyre.inventory.int("label_value", default=1)
    labelValue.meta["tip"] = "Value of label associated with initial condition patch, usually the material label value."

    db = pythia.pyre.inventory.facility("db", family="spatial_database", factory=SimpleDB)
    db.meta["tip"] = "Spatial database with values for initial condition."

    def __init__(self, name="initialconditionspatch"):
        """Constructor.
        """
        InitialCondition.__init__(self, name)

    def preinitialize(self, problem):
        """Setup initial conditions.
        """
        InitialCondition.preinitialize(self, problem)
        ModuleInitialCondition.setLabelName(self, self.labelName)
        ModuleInitialCondition.setLabelValue(self, self.labelValue)
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
    """Factory associated with InitialConditionPatch.
    """
    return InitialConditionPatch()


# End of file
