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
# @brief Python abstract basae class for specifying initial conditions.
#
# Factory: initial_conditions.

from pylith.utils.PetscComponent import PetscComponent
from .problems import InitialCondition as ModuleInitialCondition


class InitialCondition(PetscComponent, ModuleInitialCondition):
    """
    Python abstract base class for specifying initial conditions.

    INVENTORY

    Properties
      - subfields Names of solution subfields for initial condition.

    Facilities
      - None
    """

    import pythia.pyre.inventory

    subfields = pythia.pyre.inventory.list("subfields", default=["displacement"])
    subfields.meta["tip"] = "Names of solution subfields for initial condition."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="initialconditions"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="initial_conditions")
        return

    def preinitialize(self, problem):
        """
        Setup initial conditions.
        """
        self._createModuleObj()
        ModuleInitialCondition.setIdentifier(self, self.aliases[-1])
        ModuleInitialCondition.setSubfields(self, self.subfields)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        raise NotImplementedError("Please implement _createModuleOb() in derived class.")

# End of file
