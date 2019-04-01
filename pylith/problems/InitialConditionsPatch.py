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
# @file pylith/problems/InitialConditionsPatch.py
#
# @brief Python object for specifying initial conditions over a portion of the domain (patch).
#
# Factory: initial_conditions.

from pylith.problems.InitialConditions import InitialConditions
from .problems import InitialConditionsPatch as ModuleInitialConditions


def validateLabel(value):
    """
    Validate label for group/nodeset/pset.
    """
    if 0 == len(value):
        raise ValueError("Label for initial condition group/nodeset/pset in mesh not specified.")
    return value


class InitialConditionsPatch(InitialConditions, ModuleInitialConditions):
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
        InitialConditions.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Setup initial conditions.
        """
        InitialConditions.preinitialize(self, mesh)

        ModuleInitialConditions.setMarkerLabel(self, self.label)
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
    Factory associated with InitialConditionsPatch.
    """
    return InitialConditionsPatch()


# End of file
