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
# @file pylith/problems/InitialConditionPatch.py
#
# @brief Python object for specifying initial conditions over a portion of the domain (patch).
#
# Factory: initial_conditions.

from pylith.problems.InitialCondition import InitialCondition
from .problems import InitialConditionPatch as ModuleInitialCondition


def validateLabel(value):
    """
    Validate label for group/nodeset/pset.
    """
    if 0 == len(value):
        raise ValueError("Label for initial condition group/nodeset/pset in mesh not specified.")
    return value


class InitialConditionPatch(InitialCondition, ModuleInitialCondition):
    """
    Python object for specifying initial conditions over a portion of the domain (patch).

    INVENTORY

    Properties
      - *label* Label marking subdomain.

    Facilities
      - *db* Spatial database with values for initial conditions.
    """

    import pyre.inventory

    label = pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta["tip"] = "Label identifier for boundary."

    db = pyre.inventory.facility("db", default=SimpleDB, factory="spatial_database")
    db.meta["tip"] = "Spatial database with values for initial condition."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="initialconditionspatch"):
        """
        Constructor.
        """
        InitialCondition.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Setup initial conditions.
        """
        InitialCondition.preinitialize(self, mesh)

        ModuleInitialCondition.setMarkerLabel(self, self.label)
        ModuleInitialCondition.setDB(self, self.db)
        return

    def _configure(self):
        """
        Setup members using inventory.
        """
        InitialCondition._configure(self)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        ModuleInitialConditionPatch.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def initial_conditions():
    """
    Factory associated with InitialConditionPatch.
    """
    return InitialConditionPatch()


# End of file
